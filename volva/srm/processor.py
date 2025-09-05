"""
SRM (Standard Reference Material) data processing module.

Handles loading, validation, and processing of SRM diffraction patterns
for instrumental broadening correction in XRD analysis.
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, NamedTuple
import logging

from .standards import StandardReferenceMaterial, get_srm
from ..data.processor import XRDDataProcessor
from ..peaks.detector import PeakDetector
from ..core.exceptions import DataLoadError as VolvaDataError, PeakDetectionError as VolvaPeakError
from ..core.constants import WAVELENGTHS


class SRMPeakMatch(NamedTuple):
    """Matched SRM peak information."""
    observed_position: float  # 2θ in degrees
    expected_position: float  # 2θ in degrees
    hkl: Tuple[int, int, int]
    intensity: float
    fwhm: float
    position_error: float  # degrees


class SRMProcessingResult(NamedTuple):
    """Complete SRM processing results."""
    srm_name: str
    wavelength: float
    matched_peaks: List[SRMPeakMatch]
    unmatched_peaks: List[Tuple[float, float]]  # (position, intensity)
    missing_expected: List[Tuple[float, Tuple[int, int, int]]]  # (position, hkl)
    measurement_conditions: Dict[str, any]
    quality_score: float


class SRMProcessor:
    """Processes SRM diffraction data and matches with theoretical peaks."""

    def __init__(self, peak_detector: Optional[PeakDetector] = None):
        """
        Initialize SRM processor.

        Args:
            peak_detector: Peak detection instance (creates default if None)
        """
        self.peak_detector = peak_detector or PeakDetector()
        self.logger = logging.getLogger(__name__)

    def load_srm_data(self, file_path: str) -> Tuple[np.ndarray, np.ndarray, Dict[str, any]]:
        """
        Load SRM diffraction data from file.

        Args:
            file_path: Path to SRM data file

        Returns:
            Tuple of (two_theta, intensity, metadata)
        """
        try:
            processor = XRDDataProcessor()
            two_theta, intensity, message = processor.load_xrd_file(file_path)
            data = {
                'two_theta': two_theta,
                'intensity': intensity,
                'metadata': {}
            }
            # Validate SRM data quality
            self._validate_srm_data(data['two_theta'], data['intensity'])

            return data['two_theta'], data['intensity'], data.get('metadata', {})

        except Exception as e:
            raise VolvaDataError(f"Failed to load SRM data: {str(e)}")

    def _validate_srm_data(self, two_theta: np.ndarray, intensity: np.ndarray) -> None:
        """Validate SRM data quality and completeness."""
        if len(two_theta) < 100:
            raise VolvaDataError("SRM data too sparse (< 100 points)")

        if np.max(intensity) < 100:
            self.logger.warning("SRM data has low maximum intensity")

        # Check for reasonable 2θ range
        if np.max(two_theta) - np.min(two_theta) < 20:
            self.logger.warning("SRM data has limited angular range")

        # Check for negative intensities
        if np.any(intensity < 0):
            self.logger.warning("SRM data contains negative intensities")

    def identify_srm_material(self, two_theta: np.ndarray, intensity: np.ndarray,
                              wavelength: float, candidate_materials: Optional[List[str]] = None) -> str:
        """
        Automatically identify SRM material based on peak positions.

        Args:
            two_theta: 2θ positions in degrees
            intensity: Intensity values
            wavelength: X-ray wavelength in Å
            candidate_materials: List of materials to consider (all if None)

        Returns:
            Name of best matching SRM material
        """
        if candidate_materials is None:
            from .standards import list_available_srms
            candidate_materials = list_available_srms()

        # Detect peaks in the data
        try:
            peaks = self.peak_detector.detect_peaks(
                two_theta,
                intensity,
                prominence_factor=0.15,
                min_width_deg=0.05,
                use_advanced_fitting=False
            )
        except Exception as e:
            raise VolvaPeakError(f"Failed to detect peaks in SRM data: {str(e)}")

        if len(peaks) < 3:
            raise VolvaPeakError("Insufficient peaks detected for SRM identification")

        # Score each candidate material
        best_match = None
        best_score = 0

        for material_name in candidate_materials:
            srm = get_srm(material_name)
            if srm is None:
                continue

            score = self._score_srm_match(peaks, srm, wavelength)
            if score > best_score:
                best_score = score
                best_match = material_name

        if best_match is None or best_score < 0.5:
            raise VolvaPeakError("Could not identify SRM material from data")

        self.logger.info(f"Identified SRM as {best_match} (score: {best_score:.3f})")
        return best_match

    def _score_srm_match(self, observed_peaks: List, srm: StandardReferenceMaterial,
                         wavelength: float) -> float:
        """Score how well observed peaks match expected SRM peaks."""
        expected_positions = srm.get_peak_positions(wavelength, min_intensity=5.0)

        if len(expected_positions) == 0:
            return 0.0

        matches = 0
        total_expected = len(expected_positions)

        for exp_pos in expected_positions:
            # Look for observed peak within tolerance
            for peak in observed_peaks:
                peak_pos = peak.get('position', peak.get('two_theta', 0))
                if abs(peak_pos - exp_pos) <= 0.2:  # 0.2° tolerance
                    matches += 1
                    break

        return matches / total_expected

    def process_srm_data(self, two_theta: np.ndarray, intensity: np.ndarray,
                         srm_name: str, wavelength: float,
                         position_tolerance: float = 0.15,
                         min_peak_intensity: float = 5.0) -> SRMProcessingResult:
        """
        Process SRM data and match peaks with theoretical positions.

        Args:
            two_theta: 2θ positions in degrees
            intensity: Intensity values
            srm_name: Name of the SRM material
            wavelength: X-ray wavelength in Å
            position_tolerance: Peak matching tolerance in degrees
            min_peak_intensity: Minimum relative intensity for theoretical peaks

        Returns:
            Complete SRM processing results
        """
        # Get SRM reference data
        srm = get_srm(srm_name)
        if srm is None:
            raise VolvaDataError(f"Unknown SRM material: {srm_name}")

        # Detect peaks in observed data
        try:
            observed_peaks = self.peak_detector.detect_peaks(
                two_theta,
                intensity,
                prominence_factor=0.15,
                min_width_deg=0.05,
                use_advanced_fitting=False
            )
        except Exception as e:
            raise VolvaPeakError(f"Failed to detect peaks in SRM data: {str(e)}")

        # Get expected peak positions
        expected_positions = srm.get_peak_positions(wavelength, min_peak_intensity)

        # Match observed peaks with expected positions
        matched_peaks = []
        unmatched_peaks = []
        used_observed = set()

        for peak in observed_peaks:
            peak_pos = peak.get('position', peak.get('two_theta', 0))
            peak_intensity = peak.get('intensity', peak.get('height', 0))
            peak_fwhm = peak.get('fwhm', 0)

            # Find closest expected peak
            best_match = None
            min_error = float('inf')

            for exp_pos in expected_positions:
                error = abs(peak_pos - exp_pos)
                if error <= position_tolerance and error < min_error:
                    min_error = error
                    best_match = exp_pos

            if best_match is not None:
                # Get reflection info for the matched peak
                refl_info = srm.get_reflection_info(best_match, wavelength, position_tolerance)
                hkl = refl_info.hkl if refl_info else (0, 0, 0)

                matched_peaks.append(SRMPeakMatch(
                    observed_position=peak_pos,
                    expected_position=best_match,
                    hkl=hkl,
                    intensity=peak_intensity,
                    fwhm=peak_fwhm,
                    position_error=min_error
                ))
                used_observed.add(id(peak))
            else:
                unmatched_peaks.append((peak_pos, peak_intensity))

        # Find missing expected peaks
        matched_positions = {match.expected_position for match in matched_peaks}
        missing_expected = []

        for exp_pos in expected_positions:
            if exp_pos not in matched_positions:
                refl_info = srm.get_reflection_info(exp_pos, wavelength, position_tolerance)
                hkl = refl_info.hkl if refl_info else (0, 0, 0)
                missing_expected.append((exp_pos, hkl))

        # Calculate quality score
        quality_score = self._calculate_quality_score(matched_peaks, unmatched_peaks,
                                                      missing_expected, expected_positions)

        # Extract measurement conditions
        measurement_conditions = {
            'wavelength': wavelength,
            'angular_range': (np.min(two_theta), np.max(two_theta)),
            'step_size': np.mean(np.diff(two_theta)) if len(two_theta) > 1 else 0.0,
            'total_points': len(two_theta),
            'max_intensity': np.max(intensity)
        }

        return SRMProcessingResult(
            srm_name=srm_name,
            wavelength=wavelength,
            matched_peaks=matched_peaks,
            unmatched_peaks=unmatched_peaks,
            missing_expected=missing_expected,
            measurement_conditions=measurement_conditions,
            quality_score=quality_score
        )

    def _calculate_quality_score(self, matched_peaks: List[SRMPeakMatch],
                                 unmatched_peaks: List[Tuple[float, float]],
                                 missing_expected: List[Tuple[float, Tuple[int, int, int]]],
                                 expected_positions: np.ndarray) -> float:
        """Calculate overall quality score for SRM processing."""
        if len(expected_positions) == 0:
            return 0.0

        # Base score from matching rate
        match_rate = len(matched_peaks) / len(expected_positions)

        # Penalty for position errors
        if matched_peaks:
            avg_position_error = np.mean([match.position_error for match in matched_peaks])
            error_penalty = min(avg_position_error / 0.1, 1.0)  # Normalize to 0.1° reference
        else:
            error_penalty = 1.0

        # Penalty for unmatched peaks (potential contamination)
        contamination_penalty = min(len(unmatched_peaks) * 0.1, 0.5)

        # Final score
        quality_score = match_rate * (1 - error_penalty * 0.3) * (1 - contamination_penalty)
        return max(0.0, min(1.0, quality_score))

    def validate_measurement_conditions(self, srm_conditions: Dict[str, any],
                                        sample_conditions: Dict[str, any],
                                        tolerance: Dict[str, float] = None) -> Dict[str, bool]:
        """
        Validate that SRM and sample were measured under similar conditions.

        Args:
            srm_conditions: SRM measurement conditions
            sample_conditions: Sample measurement conditions
            tolerance: Tolerance values for each parameter

        Returns:
            Dictionary of validation results for each parameter
        """
        if tolerance is None:
            tolerance = {
                'wavelength': 0.0001,  # Very strict for wavelength
                'step_size': 0.001,  # 0.001° tolerance for step size
                'max_intensity': 0.2  # 20% tolerance for intensity scale
            }

        validation_results = {}

        for param in ['wavelength', 'step_size']:
            if param in srm_conditions and param in sample_conditions:
                srm_val = srm_conditions[param]
                sample_val = sample_conditions[param]
                diff = abs(srm_val - sample_val)
                validation_results[param] = diff <= tolerance.get(param, float('inf'))
            else:
                validation_results[param] = False

        # Special handling for angular range - check overlap
        if ('angular_range' in srm_conditions and
                'angular_range' in sample_conditions):
            srm_range = srm_conditions['angular_range']
            sample_range = sample_conditions['angular_range']

            # Check if ranges have reasonable overlap
            overlap_start = max(srm_range[0], sample_range[0])
            overlap_end = min(srm_range[1], sample_range[1])
            overlap_size = max(0, overlap_end - overlap_start)

            min_range_size = min(srm_range[1] - srm_range[0],
                                 sample_range[1] - sample_range[0])

            validation_results['angular_range'] = overlap_size >= 0.7 * min_range_size
        else:
            validation_results['angular_range'] = False

        return validation_results