"""
Instrumental broadening correction using SRM calibration data.

Applies SRM-derived instrumental profiles to correct sample peak widths
for accurate crystallite size and microstrain analysis.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, NamedTuple
import warnings
from pathlib import Path
import json

from .fitting import InstrumentalFitResult, InstrumentalProfileFitter
from .processor import SRMProcessingResult
from ..core.exceptions import VolvaCorrectionError
from ..peaks.models import PeakFitResult


class CorrectedPeak(NamedTuple):
    """Peak data with instrumental correction applied."""
    position: float  # 2θ in degrees
    observed_fwhm: float  # Raw FWHM in degrees
    instrumental_fwhm: float  # Instrumental FWHM at this position
    corrected_fwhm: float  # Sample-only FWHM in degrees
    intensity: float
    uncertainty: float  # Uncertainty in corrected FWHM
    correction_method: str
    deconvolution_valid: bool


class CorrectionCalibration(NamedTuple):
    """Complete instrumental calibration data."""
    srm_name: str
    wavelength: float
    fit_result: InstrumentalFitResult
    measurement_date: str
    instrument_config: Dict[str, any]
    validity_range: Tuple[float, float]  # 2θ range where calibration is valid


class InstrumentalCorrector:
    """Applies instrumental broadening corrections to sample peaks."""

    def __init__(self):
        """Initialize the corrector."""
        self.fitter = InstrumentalProfileFitter()
        self.current_calibration: Optional[CorrectionCalibration] = None

    def create_calibration(self, srm_result: SRMProcessingResult,
                           instrument_config: Dict[str, any],
                           function_name: Optional[str] = None,
                           auto_select_function: bool = True) -> CorrectionCalibration:
        """
        Create instrumental calibration from SRM data.

        Args:
            srm_result: Processed SRM measurement results
            instrument_config: Instrument configuration parameters
            function_name: Specific function to use (auto-select if None)
            auto_select_function: Whether to automatically select best function

        Returns:
            Complete calibration object
        """
        if len(srm_result.matched_peaks) < 3:
            raise VolvaCorrectionError("Insufficient SRM peaks for calibration (need ≥3)")

        # Determine best function if not specified
        if function_name is None and auto_select_function:
            try:
                fit_results = self.fitter.compare_functions(srm_result)
                function_name = self.fitter.select_best_function(fit_results, 'r_squared')
            except Exception as e:
                raise VolvaCorrectionError(f"Failed to determine best function: {str(e)}")
        elif function_name is None:
            function_name = 'Caglioti'  # Default fallback

        # Fit the instrumental profile
        try:
            fit_result = self.fitter.fit_instrumental_profile(srm_result, function_name)
        except Exception as e:
            raise VolvaCorrectionError(f"Failed to fit instrumental profile: {str(e)}")

        # Validate fit quality
        validation = self.fitter.validate_fit_quality(fit_result)
        if not validation['overall_ok']:
            warnings.warn("Instrumental fit quality is poor - corrections may be inaccurate")

        # Determine validity range
        peak_positions = [peak.observed_position for peak in srm_result.matched_peaks]
        validity_range = (min(peak_positions) - 5, max(peak_positions) + 5)

        # Create calibration
        from datetime import datetime
        calibration = CorrectionCalibration(
            srm_name=srm_result.srm_name,
            wavelength=srm_result.wavelength,
            fit_result=fit_result,
            measurement_date=datetime.now().isoformat(),
            instrument_config=instrument_config,
            validity_range=validity_range
        )

        self.current_calibration = calibration
        return calibration

    def apply_correction(self, peaks: List[PeakFitResult],
                         calibration: Optional[CorrectionCalibration] = None,
                         method: str = 'deconvolution') -> List[CorrectedPeak]:
        """
        Apply instrumental correction to sample peaks.

        Args:
            peaks: List of fitted sample peaks
            calibration: Instrumental calibration (uses current if None)
            method: Correction method ('deconvolution', 'subtraction', or 'both')

        Returns:
            List of corrected peaks
        """
        if calibration is None:
            calibration = self.current_calibration

        if calibration is None:
            raise VolvaCorrectionError("No calibration available - create one first")

        corrected_peaks = []

        for peak in peaks:
            try:
                corrected_peak = self._correct_single_peak(peak, calibration, method)
                corrected_peaks.append(corrected_peak)
            except Exception as e:
                warnings.warn(f"Failed to correct peak at {peak.position:.2f}°: {str(e)}")
                continue

        return corrected_peaks

    def _correct_single_peak(self, peak: PeakFitResult,
                             calibration: CorrectionCalibration,
                             method: str) -> CorrectedPeak:
        """Correct a single peak for instrumental broadening."""
        position = peak.position
        observed_fwhm = peak.fwhm

        # Check if position is within calibration range
        if not (calibration.validity_range[0] <= position <= calibration.validity_range[1]):
            warnings.warn(f"Peak at {position:.2f}° is outside calibration range")

        # Get instrumental FWHM at this position
        instrumental_fwhm = self.fitter.predict_instrumental_fwhm(
            np.array([position]), calibration.fit_result
        )[0]

        # Apply correction based on method
        if method == 'deconvolution':
            corrected_fwhm, uncertainty, valid = self._deconvolution_correction(
                observed_fwhm, instrumental_fwhm
            )
        elif method == 'subtraction':
            corrected_fwhm, uncertainty, valid = self._subtraction_correction(
                observed_fwhm, instrumental_fwhm
            )
        elif method == 'both':
            # Try deconvolution first, fall back to subtraction
            corrected_fwhm, uncertainty, valid = self._deconvolution_correction(
                observed_fwhm, instrumental_fwhm
            )
            if not valid:
                corrected_fwhm, uncertainty, valid = self._subtraction_correction(
                    observed_fwhm, instrumental_fwhm
                )
                method = 'subtraction'
        else:
            raise VolvaCorrectionError(f"Unknown correction method: {method}")

        return CorrectedPeak(
            position=position,
            observed_fwhm=observed_fwhm,
            instrumental_fwhm=instrumental_fwhm,
            corrected_fwhm=corrected_fwhm,
            intensity=peak.intensity,
            uncertainty=uncertainty,
            correction_method=method,
            deconvolution_valid=valid and method != 'subtraction'
        )

    def _deconvolution_correction(self, observed_fwhm: float,
                                  instrumental_fwhm: float) -> Tuple[float, float, bool]:
        """
        Deconvolve instrumental broadening assuming Gaussian profiles.

        For Gaussian profiles: FWHM_sample² = FWHM_observed² - FWHM_instrumental²
        """
        fwhm_squared_diff = observed_fwhm ** 2 - instrumental_fwhm ** 2

        if fwhm_squared_diff <= 0:
            # Invalid deconvolution - instrumental broadening dominates
            return observed_fwhm * 0.1, observed_fwhm * 0.5, False

        corrected_fwhm = np.sqrt(fwhm_squared_diff)

        # Estimate uncertainty using error propagation
        # σ(FWHM_sample) ≈ |∂FWHM_sample/∂FWHM_obs| * σ(FWHM_obs)
        # Assuming 2% uncertainty in observed FWHM
        obs_uncertainty = observed_fwhm * 0.02
        inst_uncertainty = instrumental_fwhm * 0.05  # 5% instrumental uncertainty

        derivative = observed_fwhm / corrected_fwhm
        uncertainty = np.sqrt((derivative * obs_uncertainty) ** 2 +
                              (instrumental_fwhm * inst_uncertainty / corrected_fwhm) ** 2)

        return corrected_fwhm, uncertainty, True

    def _subtraction_correction(self, observed_fwhm: float,
                                instrumental_fwhm: float) -> Tuple[float, float, bool]:
        """
        Simple subtraction correction: FWHM_sample = FWHM_observed - FWHM_instrumental

        Less accurate than deconvolution but always gives a result.
        """
        corrected_fwhm = observed_fwhm - instrumental_fwhm

        if corrected_fwhm <= 0:
            # Still invalid - use minimum reasonable value
            corrected_fwhm = observed_fwhm * 0.05
            uncertainty = observed_fwhm * 0.5
            return corrected_fwhm, uncertainty, False

        # Simple uncertainty propagation for subtraction
        uncertainty = np.sqrt((observed_fwhm * 0.02) ** 2 + (instrumental_fwhm * 0.05) ** 2)

        return corrected_fwhm, uncertainty, True

    def estimate_correction_quality(self, corrected_peaks: List[CorrectedPeak]) -> Dict[str, any]:
        """
        Assess the quality of instrumental corrections applied.

        Args:
            corrected_peaks: List of corrected peaks

        Returns:
            Dictionary with quality metrics
        """
        if not corrected_peaks:
            return {'quality': 'no_data', 'metrics': {}}

        # Count successful deconvolutions
        valid_deconvolutions = sum(1 for p in corrected_peaks if p.deconvolution_valid)
        deconvolution_rate = valid_deconvolutions / len(corrected_peaks)

        # Calculate correction magnitudes
        correction_ratios = []
        for peak in corrected_peaks:
            if peak.observed_fwhm > 0:
                ratio = peak.corrected_fwhm / peak.observed_fwhm
                correction_ratios.append(ratio)

        if correction_ratios:
            mean_ratio = np.mean(correction_ratios)
            std_ratio = np.std(correction_ratios)
        else:
            mean_ratio = std_ratio = np.nan

        # Assess relative uncertainties
        relative_uncertainties = []
        for peak in corrected_peaks:
            if peak.corrected_fwhm > 0:
                rel_unc = peak.uncertainty / peak.corrected_fwhm
                relative_uncertainties.append(rel_unc)

        mean_rel_uncertainty = np.mean(relative_uncertainties) if relative_uncertainties else np.nan

        # Overall quality assessment
        if deconvolution_rate > 0.8 and mean_ratio > 0.5 and mean_rel_uncertainty < 0.2:
            quality = 'excellent'
        elif deconvolution_rate > 0.6 and mean_ratio > 0.3 and mean_rel_uncertainty < 0.4:
            quality = 'good'
        elif deconvolution_rate > 0.3 and mean_ratio > 0.1:
            quality = 'fair'
        else:
            quality = 'poor'

        metrics = {
            'deconvolution_rate': deconvolution_rate,
            'mean_correction_ratio': mean_ratio,
            'std_correction_ratio': std_ratio,
            'mean_relative_uncertainty': mean_rel_uncertainty,
            'total_peaks': len(corrected_peaks),
            'valid_deconvolutions': valid_deconvolutions
        }

        return {'quality': quality, 'metrics': metrics}

    def save_calibration(self, calibration: CorrectionCalibration, filepath: str) -> None:
        """Save calibration to file."""
        try:
            # Convert to serializable format
            data = {
                'srm_name': calibration.srm_name,
                'wavelength': calibration.wavelength,
                'fit_result': {
                    'function_name': calibration.fit_result.function_name,
                    'parameters': calibration.fit_result.parameters,
                    'parameter_errors': calibration.fit_result.parameter_errors,
                    'r_squared': calibration.fit_result.r_squared,
                    'chi_squared': calibration.fit_result.chi_squared,
                    'two_theta_range': calibration.fit_result.two_theta_range
                },
                'measurement_date': calibration.measurement_date,
                'instrument_config': calibration.instrument_config,
                'validity_range': calibration.validity_range
            }

            with open(filepath, 'w') as f:
                json.dump(data, f, indent=2)

        except Exception as e:
            raise VolvaCorrectionError(f"Failed to save calibration: {str(e)}")

    def load_calibration(self, filepath: str) -> CorrectionCalibration:
        """Load calibration from file."""
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)

            # Reconstruct fit result
            from .fitting import InstrumentalFitResult
            fit_result = InstrumentalFitResult(
                function_name=data['fit_result']['function_name'],
                parameters=data['fit_result']['parameters'],
                parameter_errors=data['fit_result']['parameter_errors'],
                r_squared=data['fit_result']['r_squared'],
                chi_squared=data['fit_result']['chi_squared'],
                residuals=np.array([]),  # Not saved/needed for prediction
                fitted_values=np.array([]),  # Not saved/needed for prediction
                two_theta_range=tuple(data['fit_result']['two_theta_range'])
            )

            calibration = CorrectionCalibration(
                srm_name=data['srm_name'],
                wavelength=data['wavelength'],
                fit_result=fit_result,
                measurement_date=data['measurement_date'],
                instrument_config=data['instrument_config'],
                validity_range=tuple(data['validity_range'])
            )

            self.current_calibration = calibration
            return calibration

        except Exception as e:
            raise VolvaCorrectionError(f"Failed to load calibration: {str(e)}")

    def validate_calibration_for_sample(self, sample_conditions: Dict[str, any],
                                        calibration: Optional[CorrectionCalibration] = None) -> Dict[str, bool]:
        """
        Validate that current calibration is suitable for sample analysis.

        Args:
            sample_conditions: Sample measurement conditions
            calibration: Calibration to validate (uses current if None)

        Returns:
            Dictionary of validation results
        """
        if calibration is None:
            calibration = self.current_calibration

        if calibration is None:
            return {'no_calibration': False}

        validation = {}

        # Check wavelength match
        sample_wavelength = sample_conditions.get('wavelength', 0)
        wavelength_match = abs(sample_wavelength - calibration.wavelength) < 0.0001
        validation['wavelength_match'] = wavelength_match

        # Check angular range overlap
        sample_range = sample_conditions.get('angular_range', (0, 180))
        calibration_range = calibration.validity_range

        overlap_start = max(sample_range[0], calibration_range[0])
        overlap_end = min(sample_range[1], calibration_range[1])
        overlap_fraction = max(0, overlap_end - overlap_start) / (sample_range[1] - sample_range[0])

        validation['sufficient_overlap'] = overlap_fraction > 0.5

        # Check calibration age (warn if > 30 days old)
        from datetime import datetime, timedelta
        try:
            cal_date = datetime.fromisoformat(calibration.measurement_date)
            age_ok = (datetime.now() - cal_date) < timedelta(days=30)
            validation['calibration_recent'] = age_ok
        except:
            validation['calibration_recent'] = False

        # Check instrument configuration consistency
        config_match = True
        for key in ['voltage', 'current', 'detector']:
            if key in calibration.instrument_config and key in sample_conditions:
                if calibration.instrument_config[key] != sample_conditions[key]:
                    config_match = False
                    break
        validation['instrument_config_match'] = config_match

        validation['overall_valid'] = all([
            validation.get('wavelength_match', False),
            validation.get('sufficient_overlap', False),
            validation.get('instrument_config_match', False)
        ])

        return validation