"""
Volva XRD Analysis Tool - Analysis Base Classes (SRM Enhanced)
Abstract base class and common functionality for analysis methods with SRM support
"""

import numpy as np
import logging
from abc import ABC, abstractmethod
from typing import List, Dict, Optional, Tuple, Union

from ..core import AnalysisError, InsufficientDataError, ValidationError, VolvaConfig
from ..srm.correction import CorrectedPeak, CorrectionCalibration

logger = logging.getLogger(__name__)


class AnalysisMethod(ABC):
    """Abstract base class for analysis methods with SRM correction support"""

    def __init__(self, wavelength_nm: float, config: VolvaConfig):
        self.wavelength = wavelength_nm
        self.config = config
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.srm_correction_recommended = True  # Most methods benefit from SRM correction

    @abstractmethod
    def analyze(self, peaks_data: List[Dict], **kwargs) -> Dict:
        """Perform analysis on peak data"""
        pass

    def analyze_with_srm(self, peaks_data: List[Dict],
                        corrected_peaks: Optional[List[CorrectedPeak]] = None,
                        srm_calibration: Optional[CorrectionCalibration] = None,
                        **kwargs) -> Dict:
        """
        Enhanced analysis method that can use SRM-corrected peak data.

        Args:
            peaks_data: Original peak fitting results
            corrected_peaks: SRM-corrected peak data (preferred if available)
            srm_calibration: SRM calibration information
            **kwargs: Additional analysis parameters

        Returns:
            Enhanced analysis results with SRM correction tracking
        """
        # Determine if SRM correction is being used
        using_srm = corrected_peaks is not None

        if using_srm:
            # Convert corrected peaks to standard format for analysis
            enhanced_peaks_data = self._convert_corrected_peaks_to_analysis_format(
                corrected_peaks, peaks_data)
            srm_info = {
                'srm_corrected': True,
                'srm_material': srm_calibration.srm_name if srm_calibration else "Unknown",
                'correction_quality': self._assess_srm_correction_quality(corrected_peaks)
            }
        else:
            enhanced_peaks_data = peaks_data
            srm_info = {
                'srm_corrected': False,
                'srm_material': None,
                'correction_quality': None
            }

            # Add warning if SRM correction is recommended
            if self.srm_correction_recommended:
                self.logger.info("SRM correction not applied - results may have systematic errors from instrumental broadening")

        # Perform standard analysis
        result = self.analyze(enhanced_peaks_data, **kwargs)

        # Add SRM tracking information to results
        result['srm_info'] = srm_info

        # Add recommendations based on SRM usage
        if not result.get('assessment'):
            result['assessment'] = {}

        result['assessment'].update(self._generate_srm_recommendations(using_srm, srm_info))

        return result

    def _convert_corrected_peaks_to_analysis_format(self, corrected_peaks: List[CorrectedPeak],
                                                   original_peaks: List[Dict]) -> List[Dict]:
        """Convert SRM-corrected peaks back to standard analysis format."""
        enhanced_peaks = []

        for i, corr_peak in enumerate(corrected_peaks):
            # Find matching original peak by position
            original_peak = None
            for orig in original_peaks:
                if abs(orig['2theta'] - corr_peak.position) < 0.05:  # 0.05° tolerance
                    original_peak = orig.copy()
                    break

            if original_peak is None:
                # Fallback: create minimal peak data
                original_peak = {
                    '2theta': corr_peak.position,
                    'intensity': corr_peak.intensity,
                    'r_squared': 0.85
                }

            # Replace FWHM with corrected value and add SRM tracking
            enhanced_peak = original_peak.copy()
            enhanced_peak['FWHM_obs'] = corr_peak.observed_fwhm
            enhanced_peak['FWHM_corrected'] = corr_peak.corrected_fwhm  # This is the key change
            enhanced_peak['FWHM_inst'] = corr_peak.instrumental_fwhm
            enhanced_peak['srm_corrected'] = True
            enhanced_peak['srm_correction_method'] = corr_peak.correction_method
            enhanced_peak['srm_correction_valid'] = corr_peak.deconvolution_valid
            enhanced_peak['srm_correction_uncertainty'] = corr_peak.uncertainty

            enhanced_peaks.append(enhanced_peak)

        return enhanced_peaks

    def _assess_srm_correction_quality(self, corrected_peaks: List[CorrectedPeak]) -> Dict:
        """Assess the quality of SRM corrections applied."""
        if not corrected_peaks:
            return {'quality': 'no_data'}

        valid_corrections = sum(1 for p in corrected_peaks if p.deconvolution_valid)
        correction_rate = valid_corrections / len(corrected_peaks)

        # Calculate correction magnitudes
        correction_ratios = [p.corrected_fwhm / p.observed_fwhm for p in corrected_peaks
                           if p.observed_fwhm > 0]
        mean_correction = np.mean(correction_ratios) if correction_ratios else 0

        # Assess uncertainty levels
        relative_uncertainties = [p.uncertainty / p.corrected_fwhm for p in corrected_peaks
                                if p.corrected_fwhm > 0]
        mean_uncertainty = np.mean(relative_uncertainties) if relative_uncertainties else 0

        if correction_rate > 0.8 and mean_correction > 0.6 and mean_uncertainty < 0.15:
            quality = 'excellent'
        elif correction_rate > 0.6 and mean_correction > 0.4:
            quality = 'good'
        elif correction_rate > 0.3:
            quality = 'fair'
        else:
            quality = 'poor'

        return {
            'quality': quality,
            'correction_rate': correction_rate,
            'mean_correction_factor': mean_correction,
            'mean_relative_uncertainty': mean_uncertainty,
            'total_peaks': len(corrected_peaks)
        }

    def _generate_srm_recommendations(self, using_srm: bool, srm_info: Dict) -> Dict:
        """Generate recommendations based on SRM usage."""
        recommendations = []
        warnings = []

        if using_srm:
            quality = srm_info.get('correction_quality', {}).get('quality', 'unknown')
            material = srm_info.get('srm_material', 'Unknown')

            if quality == 'excellent':
                recommendations.append(f"Excellent SRM correction using {material}")
            elif quality == 'good':
                recommendations.append(f"Good SRM correction using {material}")
            elif quality == 'fair':
                warnings.append(f"Fair SRM correction quality - verify {material} calibration")
            else:
                warnings.append(f"Poor SRM correction quality - results may be unreliable")
        else:
            if self.srm_correction_recommended:
                warnings.append("SRM correction not applied - systematic instrumental broadening may affect results")
                recommendations.append("Consider measuring standard reference material (LaB6 or Si) for instrumental correction")

        return {
            'srm_recommendations': recommendations,
            'srm_warnings': warnings
        }

    def _validate_peaks(self, peaks_data: List[Dict], min_peaks: int = 1):
        """Enhanced peak data validation with SRM awareness"""
        if not peaks_data:
            raise InsufficientDataError('No peaks provided')
        if len(peaks_data) < min_peaks:
            raise InsufficientDataError(f'Need at least {min_peaks} peaks, got {len(peaks_data)}')

        # Determine which FWHM field to use
        fwhm_field = 'FWHM_corrected' if peaks_data[0].get('srm_corrected', False) else 'FWHM_obs'

        required_fields = ['2theta', fwhm_field]
        for i, peak in enumerate(peaks_data):
            for field in required_fields:
                if field not in peak:
                    raise AnalysisError(f'Peak {i + 1} missing required field: {field}')

            # Scientific validation
            two_theta = peak['2theta']
            fwhm = peak[fwhm_field]

            # Check reasonable 2theta range
            if not (5 <= two_theta <= 180):
                self.logger.warning(f"Peak {i + 1} at unusual 2θ = {two_theta:.2f}°")

            # Check reasonable FWHM range
            if not (0.005 <= fwhm <= 5.0):  # Lowered minimum for SRM-corrected peaks
                raise ValidationError(f"Peak {i + 1} has unrealistic FWHM = {fwhm:.4f}°")

            # Check peak quality if available
            r_squared = peak.get('r_squared', 1.0)
            if r_squared < 0.5:
                self.logger.warning(f"Peak {i + 1} has poor fit quality (R² = {r_squared:.3f})")

            # SRM-specific validations
            if peak.get('srm_corrected', False):
                if not peak.get('srm_correction_valid', True):
                    self.logger.warning(f"Peak {i + 1} has invalid SRM deconvolution")

                uncertainty = peak.get('srm_correction_uncertainty', 0)
                if uncertainty > fwhm * 0.5:
                    self.logger.warning(f"Peak {i + 1} has large SRM correction uncertainty")

    def _correct_instrumental_broadening(self, fwhm_obs: float,
                                         fwhm_inst: float = None) -> float:
        """
        Correct for instrumental broadening - now SRM-aware.

        Note: This method is kept for backward compatibility, but when SRM
        corrections are available, the corrected FWHM should be used directly.
        """
        if fwhm_inst is None:
            fwhm_inst = self.config.default_instrumental_broadening

        # Warren's formula: β_sample² = β_observed² - β_instrument²
        beta_corrected_sq = max(fwhm_obs ** 2 - fwhm_inst ** 2, (fwhm_obs * 0.1) ** 2)
        beta_corrected = np.sqrt(beta_corrected_sq)

        # Log significant corrections
        correction_factor = beta_corrected / fwhm_obs
        if correction_factor < 0.8:
            self.logger.warning(f"Large instrumental correction: "
                                f"FWHM {fwhm_obs:.4f}° → {beta_corrected:.4f}° "
                                f"({correction_factor:.1%} retained)")

        return beta_corrected

    def _get_effective_fwhm(self, peak: Dict) -> float:
        """
        Get the appropriate FWHM value for analysis.

        Prioritizes SRM-corrected FWHM over instrumental correction.
        """
        if peak.get('srm_corrected', False):
            # Use SRM-corrected FWHM directly
            return peak['FWHM_corrected']
        else:
            # Fall back to traditional instrumental correction
            fwhm_obs = peak['FWHM_obs']
            fwhm_inst = peak.get('FWHM_inst', self.config.default_instrumental_broadening)
            return self._correct_instrumental_broadening(fwhm_obs, fwhm_inst)

    def _validate_crystallite_size(self, size_nm: float, method_name: str) -> bool:
        """Validate calculated crystallite size for physical reasonableness"""
        if size_nm <= 0:
            self.logger.error(f"{method_name}: Negative crystallite size")
            return False

        if size_nm < 0.5:  # Sub-nanometer
            self.logger.warning(f"{method_name}: Very small crystallite size ({size_nm:.2f} nm)")
        elif size_nm > 1000:  # Micron scale
            self.logger.warning(f"{method_name}: Very large crystallite size ({size_nm:.0f} nm)")

        if size_nm > 10000:  # Unrealistic for XRD
            self.logger.error(f"{method_name}: Unrealistic crystallite size ({size_nm:.0f} nm)")
            return False

        return True

    def _check_scherrer_limits(self, two_theta: float, fwhm_deg: float) -> Dict:
        """Check if Scherrer equation limitations are met"""
        warnings = []

        # Convert to radians for calculation
        theta_rad = np.deg2rad(two_theta / 2)
        fwhm_rad = np.deg2rad(fwhm_deg)

        # Check if peak is broad enough for reliable size determination
        if fwhm_deg < 0.03:  # Lowered threshold for SRM-corrected peaks
            warnings.append("Peak may be too sharp for reliable size analysis")

        # Check angular range - Scherrer works best at medium angles
        if two_theta < 20:
            warnings.append("Low angle peak - consider absorption corrections")
        elif two_theta > 120:
            warnings.append("High angle peak - check for systematic errors")

        # Estimate size to check if it's in the valid range for Scherrer
        approx_size_nm = (0.9 * self.wavelength * 1e9) / (fwhm_rad * np.cos(theta_rad))

        if approx_size_nm > 200:
            warnings.append(f"Estimated size ({approx_size_nm:.0f} nm) approaching Scherrer limit")

        return {
            'valid': len(warnings) == 0,
            'warnings': warnings,
            'estimated_size_nm': approx_size_nm
        }


class StrainAnalysisMethod(AnalysisMethod):
    """Base class for strain-based analysis methods with SRM support"""

    def __init__(self, wavelength_nm: float, config: VolvaConfig):
        super().__init__(wavelength_nm, config)
        self.srm_correction_recommended = True  # Strain analysis especially benefits from SRM

    def _validate_strain_analysis(self, peaks_data: List[Dict]) -> Dict:
        """Enhanced validation for strain analysis methods with SRM awareness"""
        validation_results = {
            'valid': True,
            'warnings': [],
            'recommendations': []
        }

        # Check minimum number of peaks
        if len(peaks_data) < 3:
            validation_results['valid'] = False
            validation_results['warnings'].append("Need at least 3 peaks for strain analysis")
            return validation_results

        if len(peaks_data) < 5:
            validation_results['warnings'].append("Recommend 5+ peaks for reliable strain analysis")

        # Check SRM correction status
        srm_corrected_count = sum(1 for peak in peaks_data if peak.get('srm_corrected', False))
        if srm_corrected_count > 0:
            if srm_corrected_count == len(peaks_data):
                validation_results['recommendations'].append("All peaks have SRM correction - excellent for strain analysis")
            else:
                validation_results['warnings'].append(f"Mixed SRM correction status ({srm_corrected_count}/{len(peaks_data)} corrected)")
        else:
            validation_results['warnings'].append("No SRM correction applied - strain results may be less accurate")

        # Check angular range coverage
        angles = [p['2theta'] for p in peaks_data]
        angle_range = max(angles) - min(angles)

        if angle_range < 30:
            validation_results['warnings'].append(f"Limited angular range ({angle_range:.1f}°)")
            validation_results['recommendations'].append("Use peaks across wider 2θ range")

        # Use appropriate FWHM field for analysis
        fwhm_values = [self._get_effective_fwhm(peak) for peak in peaks_data]

        # Simple check: FWHM should generally increase with angle for strain
        if len(set(np.round(fwhm_values, 4))) == 1:  # All FWHM values identical
            validation_results['warnings'].append("Constant peak widths - check for size-only broadening")

        # Check peak quality distribution
        qualities = [p.get('r_squared', 0.5) for p in peaks_data]
        avg_quality = np.mean(qualities)

        if avg_quality < 0.7:
            validation_results['warnings'].append(f"Average peak quality low (R² = {avg_quality:.3f})")

        # Check for potential texture effects
        intensities = [p['intensity'] for p in peaks_data if 'intensity' in p]
        if intensities:
            intensity_cv = np.std(intensities) / np.mean(intensities)
            if intensity_cv > 2.0:
                validation_results['warnings'].append("Large intensity variations - check for preferred orientation")
                validation_results['recommendations'].append("Consider texture corrections")

        return validation_results

    def _assess_isotropy_assumption(self, x_data: np.ndarray, y_data: np.ndarray,
                                    r_squared: float) -> Dict:
        """Assess whether isotropic strain model is appropriate"""
        assessment = {
            'likely_isotropic': True,
            'confidence': 'high',
            'warnings': [],
            'recommendations': []
        }

        # R² assessment
        if r_squared > 0.95:
            assessment['confidence'] = 'very_high'
        elif r_squared > 0.9:
            assessment['confidence'] = 'high'
        elif r_squared > 0.8:
            assessment['confidence'] = 'moderate'
            assessment['warnings'].append(f"Moderate linear fit (R² = {r_squared:.3f})")
        else:
            assessment['likely_isotropic'] = False
            assessment['confidence'] = 'low'
            assessment['warnings'].append(f"Poor linear fit (R² = {r_squared:.3f}) suggests anisotropic strain")
            assessment['recommendations'].append("Consider anisotropic strain models")

        # Check for systematic deviations
        if len(x_data) >= 4:
            # Fit line and check residual patterns
            coeffs = np.polyfit(x_data, y_data, 1)
            y_pred = np.polyval(coeffs, x_data)
            residuals = y_data - y_pred

            # Check for systematic curvature
            if len(x_data) >= 6:
                try:
                    quad_coeffs = np.polyfit(x_data, y_data, 2)
                    quad_pred = np.polyval(quad_coeffs, x_data)
                    quad_r2 = 1 - np.sum((y_data - quad_pred) ** 2) / np.sum((y_data - np.mean(y_data)) ** 2)

                    if quad_r2 > r_squared + 0.05:
                        assessment['warnings'].append("Quadratic fit significantly better - consider H-W method")
                        assessment['recommendations'].append("Try Halder-Wagner analysis")
                except:
                    pass

            # Check for outliers
            residual_std = np.std(residuals)
            outliers = np.abs(residuals) > 2.5 * residual_std
            if np.any(outliers):
                n_outliers = np.sum(outliers)
                assessment['warnings'].append(f"{n_outliers} potential outlier(s) detected")
                assessment['recommendations'].append("Check peak quality and consider data filtering")

        return assessment

    def _calculate_strain_uncertainty(self, slope: float, slope_error: float,
                                      method: str) -> Tuple[float, str]:
        """Calculate microstrain uncertainty with SRM correction awareness"""
        try:
            if not np.isfinite(slope_error) or slope_error <= 0:
                return 0.0, "uncertainty_unavailable"

            if method.lower() == 'williamson_hall':
                strain_error = (slope_error / 4) * 100  # Convert to percentage
            elif method.lower() == 'halder_wagner':
                if slope <= 0:
                    return 0.0, "negative_slope"
                a_param = getattr(self.config, 'default_hw_a', 2.0)
                strain_error = (slope_error / (2 * np.sqrt(slope * a_param))) * 100
            else:
                return 0.0, "unknown_method"

            return strain_error, "calculated"

        except Exception as e:
            logger.warning(f"Strain uncertainty calculation failed: {e}")
            return 0.0, "calculation_error"