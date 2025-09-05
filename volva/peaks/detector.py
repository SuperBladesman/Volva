"""
Volva XRD Analysis Tool - Peak Detection Module
Enhanced peak detection and fitting with multiple methods
"""

import numpy as np
from scipy.signal import find_peaks
from typing import List, Dict, Optional
import logging

from ..core import PeakDetectionError, InsufficientDataError, VolvaConfig
from ..data import XRDDataProcessor
from .models import PeakFitter

logger = logging.getLogger(__name__)


class PeakDetector:
    """Enhanced peak detection and fitting with multiple methods"""

    def __init__(self, config: VolvaConfig):
        self.config = config
        self.fitter = PeakFitter(use_lmfit=config.use_lmfit)

    def detect_peaks(self, two_theta: np.ndarray, intensity: np.ndarray,
                     prominence_factor: float = None, min_width_deg: float = None,
                     use_advanced_fitting: bool = None) -> List[Dict]:
        """Enhanced peak detection with configurable parameters"""

        # Use config defaults if not specified
        prominence_factor = prominence_factor or self.config.default_prominence_factor
        min_width_deg = min_width_deg or self.config.default_min_width_deg
        use_advanced_fitting = use_advanced_fitting if use_advanced_fitting is not None else self.config.use_lmfit

        try:
            if len(two_theta) < 10:
                raise InsufficientDataError("Need at least 10 data points for peak detection")

            # Smooth data
            smoothed = XRDDataProcessor.smooth_data(intensity)

            # Calculate detection parameters
            step_size = np.median(np.diff(two_theta))
            min_width_points = max(3, int(min_width_deg / step_size))
            min_distance = max(5, int(0.1 / step_size))

            # Adaptive threshold calculation
            baseline = np.percentile(smoothed, 10)
            noise_level = np.std(smoothed[smoothed < np.percentile(smoothed, 30)])
            peak_height = baseline + max(3 * noise_level, prominence_factor * (np.max(smoothed) - baseline))

            # Find peaks
            peaks, properties = find_peaks(
                smoothed,
                height=peak_height,
                distance=min_distance,
                width=min_width_points,
                prominence=prominence_factor * np.max(smoothed)
            )

            if len(peaks) == 0:
                logger.warning("No peaks detected")
                return []

            peaks_info = []
            for i, peak_idx in enumerate(peaks):
                if peak_idx >= len(two_theta):
                    continue

                peak_2theta = two_theta[peak_idx]
                peak_intensity = intensity[peak_idx]

                # Enhanced FWHM estimation
                if use_advanced_fitting:
                    try:
                        fwhm_estimate, r_squared, fit_info = self._fit_peak_advanced(
                            two_theta, intensity, peak_2theta)
                    except Exception as e:
                        logger.warning(f"Advanced fitting failed for peak {i + 1}: {e}")
                        fwhm_estimate, r_squared, fit_info = self._estimate_fwhm_simple(
                            two_theta, intensity, peak_2theta)
                else:
                    fwhm_estimate, r_squared, fit_info = self._estimate_fwhm_simple(
                        two_theta, intensity, peak_2theta)

                # Estimate instrumental broadening based on 2theta
                inst_broadening = self._estimate_instrumental_broadening(peak_2theta)

                # Get additional peak properties from scipy
                peak_properties = self._extract_peak_properties(properties, i)

                peaks_info.append({
                    'index': peak_idx,
                    '2theta': peak_2theta,
                    'intensity': peak_intensity,
                    'FWHM_obs': fwhm_estimate,
                    'FWHM_inst': inst_broadening,
                    'r_squared': r_squared,
                    'fit_info': fit_info,
                    'prominence': peak_properties.get('prominence', 0),
                    'width_points': peak_properties.get('width', 0),
                    'left_base': peak_properties.get('left_base', peak_idx),
                    'right_base': peak_properties.get('right_base', peak_idx)
                })

            logger.info(f"Peak detection completed: {len(peaks_info)} peaks found")
            return peaks_info

        except Exception as e:
            if isinstance(e, (InsufficientDataError, PeakDetectionError)):
                raise
            raise PeakDetectionError(f"Peak detection failed: {str(e)}")

    def _fit_peak_advanced(self, two_theta: np.ndarray, intensity: np.ndarray,
                           peak_center: float) -> tuple:
        """Advanced peak fitting using lmfit"""
        try:
            # Define fitting region
            mask = np.abs(two_theta - peak_center) <= self.config.default_fitting_window_deg
            x_region = two_theta[mask]
            y_region = intensity[mask]

            if len(x_region) < 5:
                raise ValueError("Insufficient data points in fitting region")

            return self.fitter.fit_peak_advanced(
                x_region, y_region, peak_center, self.config.default_peak_model)

        except Exception as e:
            logger.debug(f"Advanced peak fitting failed: {e}")
            return self._estimate_fwhm_simple(two_theta, intensity, peak_center)

    def _estimate_fwhm_simple(self, two_theta: np.ndarray, intensity: np.ndarray,
                              peak_center: float) -> tuple:
        """Simple FWHM estimation using half-maximum method"""
        try:
            # Find peak region
            mask = np.abs(two_theta - peak_center) <= 1.0
            x_region = two_theta[mask]
            y_region = intensity[mask]

            if len(x_region) < 5:
                return 0.1, 0.5, {'method': 'failed', 'reason': 'insufficient_data'}

            return self.fitter.fit_peak_simple(x_region, y_region, peak_center)

        except Exception:
            return 0.1, 0.3, {'method': 'error', 'reason': 'exception'}

    def _estimate_instrumental_broadening(self, two_theta: float) -> float:
        """Estimate instrumental broadening based on 2theta (Caglioti function approximation)"""
        # Simple approximation: instrumental broadening increases with angle
        base_broadening = self.config.default_instrumental_broadening
        # U*tan²θ + V*tanθ + W approximation (simplified)
        theta_rad = np.deg2rad(two_theta / 2)
        caglioti_factor = 1 + 0.1 * np.tan(theta_rad) ** 2
        return base_broadening * caglioti_factor

    def _extract_peak_properties(self, properties: Dict, peak_index: int) -> Dict:
        """Extract additional properties from scipy find_peaks result"""
        try:
            peak_props = {}

            if 'prominences' in properties and peak_index < len(properties['prominences']):
                peak_props['prominence'] = properties['prominences'][peak_index]

            if 'widths' in properties and peak_index < len(properties['widths']):
                peak_props['width'] = properties['widths'][peak_index]

            if 'left_bases' in properties and peak_index < len(properties['left_bases']):
                peak_props['left_base'] = int(properties['left_bases'][peak_index])

            if 'right_bases' in properties and peak_index < len(properties['right_bases']):
                peak_props['right_base'] = int(properties['right_bases'][peak_index])

            if 'peak_heights' in properties and peak_index < len(properties['peak_heights']):
                peak_props['peak_height'] = properties['peak_heights'][peak_index]

            return peak_props

        except Exception as e:
            logger.debug(f"Failed to extract peak properties: {e}")
            return {}

    def refine_peak_positions(self, peaks_info: List[Dict], two_theta: np.ndarray,
                              intensity: np.ndarray) -> List[Dict]:
        """Refine peak positions using more sophisticated fitting"""
        if not self.config.use_lmfit:
            return peaks_info

        try:
            refined_peaks = []

            for peak in peaks_info:
                peak_center = peak['2theta']

                # Define larger fitting region for refinement
                mask = np.abs(two_theta - peak_center) <= self.config.default_fitting_window_deg * 1.5
                x_region = two_theta[mask]
                y_region = intensity[mask]

                if len(x_region) < 10:
                    refined_peaks.append(peak)
                    continue

                try:
                    fwhm, r_squared, fit_info = self.fitter.fit_peak_advanced(
                        x_region, y_region, peak_center, self.config.default_peak_model)

                    # Update peak information with refined values
                    refined_peak = peak.copy()
                    refined_peak['FWHM_obs'] = fwhm
                    refined_peak['r_squared'] = r_squared
                    refined_peak['fit_info'] = fit_info

                    # Update center position if fitting provides it
                    if 'fit_params' in fit_info and 'center' in fit_info['fit_params']:
                        new_center = fit_info['fit_params']['center']
                        if abs(new_center - peak_center) < 0.5:  # Reasonable refinement
                            refined_peak['2theta'] = new_center
                            # Find corresponding intensity
                            center_idx = np.argmin(np.abs(two_theta - new_center))
                            refined_peak['intensity'] = intensity[center_idx]
                            refined_peak['index'] = center_idx

                    refined_peaks.append(refined_peak)

                except Exception as e:
                    logger.debug(f"Peak refinement failed for peak at {peak_center:.2f}°: {e}")
                    refined_peaks.append(peak)

            logger.info(f"Peak position refinement completed for {len(refined_peaks)} peaks")
            return refined_peaks

        except Exception as e:
            logger.warning(f"Peak refinement failed: {e}")
            return peaks_info

    def filter_peaks_by_quality(self, peaks_info: List[Dict],
                                min_r_squared: float = None) -> List[Dict]:
        """Filter peaks based on quality criteria"""
        min_r_squared = min_r_squared or self.config.min_r_squared

        try:
            filtered_peaks = []

            for peak in peaks_info:
                r_squared = peak.get('r_squared', 0.0)
                fwhm = peak.get('FWHM_obs', 0.0)
                intensity = peak.get('intensity', 0.0)

                # Quality criteria
                quality_checks = [
                    r_squared >= min_r_squared,
                    0.01 <= fwhm <= 2.0,  # Reasonable FWHM range
                    intensity > 0,
                    peak.get('prominence', 0) > 0
                ]

                if all(quality_checks):
                    filtered_peaks.append(peak)
                else:
                    logger.debug(f"Peak at {peak['2theta']:.2f}° filtered out: "
                                 f"R²={r_squared:.3f}, FWHM={fwhm:.3f}°")

            logger.info(f"Quality filtering: {len(filtered_peaks)}/{len(peaks_info)} peaks retained")
            return filtered_peaks

        except Exception as e:
            logger.warning(f"Peak quality filtering failed: {e}")
            return peaks_info

    def group_overlapping_peaks(self, peaks_info: List[Dict]) -> List[List[Dict]]:
        """Group peaks that might be overlapping"""
        if len(peaks_info) < 2:
            return [[peak] for peak in peaks_info]

        try:
            # Sort peaks by 2theta
            sorted_peaks = sorted(peaks_info, key=lambda x: x['2theta'])
            groups = []
            current_group = [sorted_peaks[0]]

            for i in range(1, len(sorted_peaks)):
                current_peak = sorted_peaks[i]
                last_peak = current_group[-1]

                # Check if peaks are close enough to be overlapping
                separation = abs(current_peak['2theta'] - last_peak['2theta'])
                avg_fwhm = (current_peak['FWHM_obs'] + last_peak['FWHM_obs']) / 2

                # If separation is less than sum of half-widths, consider overlapping
                if separation < avg_fwhm:
                    current_group.append(current_peak)
                else:
                    groups.append(current_group)
                    current_group = [current_peak]

            groups.append(current_group)

            # Log overlapping groups
            overlapping_groups = [g for g in groups if len(g) > 1]
            if overlapping_groups:
                logger.info(f"Found {len(overlapping_groups)} groups of overlapping peaks")
                for i, group in enumerate(overlapping_groups):
                    positions = [f"{p['2theta']:.2f}°" for p in group]
                    logger.debug(f"Overlap group {i + 1}: {', '.join(positions)}")

            return groups

        except Exception as e:
            logger.warning(f"Peak grouping failed: {e}")
            return [[peak] for peak in peaks_info]

    def get_peak_statistics(self, peaks_info: List[Dict]) -> Dict:
        """Calculate statistics for detected peaks"""
        if not peaks_info:
            return {}

        try:
            # Extract peak properties
            positions = [p['2theta'] for p in peaks_info]
            intensities = [p['intensity'] for p in peaks_info]
            fwhms = [p['FWHM_obs'] for p in peaks_info]
            r_squareds = [p.get('r_squared', 0.5) for p in peaks_info]
            prominences = [p.get('prominence', 0) for p in peaks_info]

            statistics = {
                'num_peaks': len(peaks_info),
                'position_range': (min(positions), max(positions)),
                'avg_intensity': np.mean(intensities),
                'std_intensity': np.std(intensities),
                'avg_fwhm': np.mean(fwhms),
                'std_fwhm': np.std(fwhms),
                'avg_r_squared': np.mean(r_squareds),
                'min_r_squared': min(r_squareds),
                'max_r_squared': max(r_squareds),
                'high_quality_peaks': sum(1 for r2 in r_squareds if r2 > 0.8),
                'medium_quality_peaks': sum(1 for r2 in r_squareds if 0.6 < r2 <= 0.8),
                'low_quality_peaks': sum(1 for r2 in r_squareds if r2 <= 0.6),
                'avg_prominence': np.mean(prominences) if prominences else 0,
                'total_intensity': sum(intensities)
            }

            return statistics

        except Exception as e:
            logger.warning(f"Peak statistics calculation failed: {e}")
            return {'num_peaks': len(peaks_info)}