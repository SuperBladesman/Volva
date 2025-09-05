"""
Volva XRD Analysis Tool - Monshi-Scherrer Analysis
Enhanced Monshi-Scherrer analysis with proper Kα doublet correction
"""

import numpy as np
from typing import List, Dict, Union, Tuple
import logging

from ..core import AnalysisError, Constants, VolvaConfig
from .base import AnalysisMethod

logger = logging.getLogger(__name__)


class MonshiScherrerAnalyzer(AnalysisMethod):
    """Enhanced Monshi-Scherrer analysis with proper Kα2 doublet handling"""

    def __init__(self, wavelength_info: Dict, config: VolvaConfig):
        # Extract primary wavelength for parent class
        ka1_wavelength = wavelength_info.get('ka1', 0.15406)
        super().__init__(ka1_wavelength, config)

        # Store complete wavelength information
        self.wavelength_info = wavelength_info
        self.ka1 = wavelength_info.get('ka1', 0.15406)
        self.ka2 = wavelength_info.get('ka2', 0.154439)
        self.ka2_ka1_ratio = wavelength_info.get('ratio', 0.5)
        self.alpha = wavelength_info.get('alpha', 2.2e-3)

        # Validate wavelength data
        self._validate_wavelength_data()

        self.logger.info(f"Monshi-Scherrer initialized: Kα1={self.ka1:.5f}nm, "
                         f"Kα2={self.ka2:.5f}nm, α={self.alpha:.3e}")

    def _validate_wavelength_data(self):
        """Validate wavelength information for physical consistency"""
        if self.ka2 <= self.ka1:
            self.logger.warning(f"Unusual wavelength order: Kα2 ({self.ka2:.5f}) <= Kα1 ({self.ka1:.5f})")

        if not (0.1 <= self.ka2_ka1_ratio <= 1.0):
            self.logger.warning(f"Unusual Kα2/Kα1 intensity ratio: {self.ka2_ka1_ratio:.3f}")

        if not (1e-4 <= self.alpha <= 1e-2):
            self.logger.warning(f"Unusual alpha parameter: {self.alpha:.3e}")

    def analyze(self, peaks_data: List[Dict], k_constant: Union[float, str] = None,
                doublet_method: str = 'correction') -> Dict:
        """
        Enhanced Monshi-Scherrer analysis with proper Kα doublet handling

        Parameters:
        - doublet_method: 'correction' (default) applies Monshi correction,
                         'deconvolution' attempts to separate Kα1/Kα2 components
        """
        self._validate_peaks(peaks_data)

        # Handle K constant
        k_value, k_description = self._process_k_constant(k_constant)

        # Validate doublet method
        if doublet_method not in ['correction', 'deconvolution']:
            self.logger.warning(f"Unknown doublet method '{doublet_method}', using 'correction'")
            doublet_method = 'correction'

        results = []
        failed_peaks = []
        validation_warnings = []
        doublet_warnings = []

        for i, peak in enumerate(peaks_data):
            try:
                peak_result = self._analyze_single_peak(peak, k_value, i + 1, doublet_method)
                if peak_result['valid']:
                    results.append(peak_result['data'])
                    validation_warnings.extend(peak_result['warnings'])
                    doublet_warnings.extend(peak_result.get('doublet_warnings', []))
                else:
                    failed_peaks.append({
                        'peak_number': i + 1,
                        'reason': peak_result['reason'],
                        '2theta': peak['2theta']
                    })

            except Exception as e:
                self.logger.warning(f"Monshi-Scherrer calculation failed for peak {i + 1}: {e}")
                failed_peaks.append({
                    'peak_number': i + 1,
                    'reason': f'Calculation error: {str(e)}',
                    '2theta': peak['2theta']
                })
                continue

        if not results:
            raise AnalysisError('No valid peaks for Monshi-Scherrer analysis')

        # Calculate statistics
        sizes = [r['size_nm'] for r in results]
        quality_weights = [r['r_squared'] for r in results]
        weighted_avg = np.average(sizes, weights=quality_weights)

        analysis_result = {
            'success': True,
            'method': 'Monshi-Scherrer',
            'doublet_method': doublet_method,
            'results': results,
            'average_size': np.mean(sizes),
            'weighted_average_size': weighted_avg,
            'std_size': np.std(sizes),
            'median_size': np.median(sizes),
            'size_range': (min(sizes), max(sizes)),
            'num_peaks': len(results),
            'failed_peaks': failed_peaks,
            'k_constant': k_value,
            'k_description': k_description,
            'ka1_wavelength': self.ka1,
            'ka2_wavelength': self.ka2,
            'ka2_ka1_ratio': self.ka2_ka1_ratio,
            'alpha_parameter': self.alpha,
            'validation_warnings': validation_warnings,
            'doublet_warnings': doublet_warnings,
            'quality_metrics': self._calculate_quality_metrics(results)
        }

        # Overall assessment
        analysis_result['assessment'] = self._assess_analysis_quality(analysis_result)

        return analysis_result

    def _process_k_constant(self, k_constant: Union[float, str, None]) -> Tuple[str, str]:
        """Process K constant with Monshi-Scherrer specific validation"""
        if k_constant is None:
            k_constant = self.config.default_scherrer_k

        if isinstance(k_constant, str):
            if k_constant not in Constants.SCHERRER_CONSTANTS:
                available = list(Constants.SCHERRER_CONSTANTS.keys())
                raise ValueError(f"Unknown Scherrer constant: '{k_constant}'. Available: {available}")
            k_data = Constants.SCHERRER_CONSTANTS[k_constant]
            k_value = k_data["value"]
            k_description = f"{k_constant} ({k_data['description']})"
        else:
            k_value = float(k_constant)
            k_description = f"Custom value: {k_value}"

        return k_value, k_description

    def _analyze_single_peak(self, peak: Dict, k_value: float, peak_number: int,
                             doublet_method: str) -> Dict:
        """Analyze single peak with proper Kα doublet treatment"""
        two_theta = peak['2theta']
        fwhm_obs = peak['FWHM_obs']
        fwhm_inst = peak.get('FWHM_inst', self.config.default_instrumental_broadening)

        warnings = []
        doublet_warnings = []

        # Check Scherrer limits
        scherrer_check = self._check_scherrer_limits(two_theta, fwhm_obs)
        warnings.extend(scherrer_check['warnings'])

        # Step 1: Correct for instrumental broadening (standard Warren correction)
        beta_corrected = self._get_effective_fwhm(peak)
        correction_factor = beta_corrected / fwhm_obs

        # Step 2: Apply doublet correction
        if doublet_method == 'correction':
            beta_monshi, doublet_info = self._apply_monshi_correction(
                beta_corrected, two_theta, peak_number)
        elif doublet_method == 'deconvolution':
            beta_monshi, doublet_info = self._apply_doublet_deconvolution(
                beta_corrected, two_theta, peak_number)
        else:
            beta_monshi = beta_corrected
            doublet_info = {'method': 'none', 'correction_applied': False}

        doublet_warnings.extend(doublet_info.get('warnings', []))

        # Step 3: Calculate crystallite size using Kα1 wavelength
        try:
            theta_rad = np.deg2rad(two_theta / 2.0)
            cos_theta = np.cos(theta_rad)

            if cos_theta <= 0:
                return {
                    'valid': False,
                    'reason': f'Invalid angle: cos(θ) = {cos_theta:.6f}',
                    'warnings': warnings,
                    'doublet_warnings': doublet_warnings
                }

            beta_monshi_rad = np.deg2rad(beta_monshi)
            crystallite_size = k_value * self.ka1 / (beta_monshi_rad * cos_theta)
            crystallite_size_nm = crystallite_size * 1e9

        except (ZeroDivisionError, ValueError) as e:
            return {
                'valid': False,
                'reason': f'Calculation error: {str(e)}',
                'warnings': warnings,
                'doublet_warnings': doublet_warnings
            }

        # Validate result
        if not self._validate_crystallite_size(crystallite_size_nm, f"Peak {peak_number}"):
            return {
                'valid': False,
                'reason': f'Unrealistic size: {crystallite_size_nm:.1f} nm',
                'warnings': warnings,
                'doublet_warnings': doublet_warnings
            }

        # Calculate uncertainty
        size_uncertainty = self._estimate_size_uncertainty(
            crystallite_size_nm, beta_monshi, beta_corrected, peak, doublet_info)

        peak_result = {
            'peak_number': peak_number,
            'two_theta': two_theta,
            'size_nm': crystallite_size_nm,
            'size_uncertainty_nm': size_uncertainty,
            'fwhm_obs': fwhm_obs,
            'fwhm_corrected': beta_corrected,
            'fwhm_monshi': beta_monshi,
            'fwhm_inst': fwhm_inst,
            'correction_factor': correction_factor,
            'doublet_correction_factor': beta_monshi / beta_corrected,
            'r_squared': peak.get('r_squared', 0.85),
            'doublet_info': doublet_info
        }

        return {
            'valid': True,
            'data': peak_result,
            'warnings': warnings,
            'doublet_warnings': doublet_warnings
        }

    def _apply_monshi_correction(self, beta_corrected: float, two_theta: float,
                                 peak_number: int) -> Tuple[float, Dict]:
        """
        Apply Monshi correction for Kα doublet broadening

        The Monshi method corrects for the fact that observed peaks are convolutions
        of Kα1 and Kα2 components. The correction removes the doublet broadening.
        """
        try:
            theta_rad = np.deg2rad(two_theta / 2.0)

            # Monshi correction: β_corrected = β_observed - α*tan(θ)
            # where α is the doublet parameter
            doublet_broadening_deg = self.alpha * np.tan(theta_rad)

            # The correction should be subtracted from the FWHM
            beta_monshi = beta_corrected - doublet_broadening_deg

            # Ensure we don't get negative or unreasonably small values
            min_beta = beta_corrected * 0.1  # Minimum 10% of corrected FWHM
            if beta_monshi < min_beta:
                warnings = [f"Large Monshi correction for peak {peak_number}: "
                            f"{doublet_broadening_deg:.4f}° removed from {beta_corrected:.4f}°"]
                beta_monshi = min_beta
            else:
                warnings = []

            # Calculate correction statistics
            correction_magnitude = doublet_broadening_deg / beta_corrected

            doublet_info = {
                'method': 'monshi_correction',
                'correction_applied': True,
                'doublet_broadening_deg': doublet_broadening_deg,
                'correction_magnitude': correction_magnitude,
                'alpha_used': self.alpha,
                'theta_rad': theta_rad,
                'warnings': warnings
            }

            # Add warnings for large corrections
            if correction_magnitude > 0.3:  # >30% correction
                doublet_info['warnings'].append(
                    f"Large doublet correction ({correction_magnitude:.1%}) - verify alpha parameter")

            return beta_monshi, doublet_info

        except Exception as e:
            self.logger.warning(f"Monshi correction failed for peak {peak_number}: {e}")
            return beta_corrected, {
                'method': 'monshi_correction',
                'correction_applied': False,
                'error': str(e),
                'warnings': ['Monshi correction failed, using instrumental-corrected FWHM']
            }

    def _apply_doublet_deconvolution(self, beta_corrected: float, two_theta: float,
                                     peak_number: int) -> Tuple[float, Dict]:
        """
        Alternative approach: Attempt to deconvolve Kα1 and Kα2 components

        This method estimates the Kα1 component width by modeling the doublet
        as two overlapping peaks with known intensity ratio and separation.
        """
        try:
            theta_rad = np.deg2rad(two_theta / 2.0)

            # Calculate Kα1-Kα2 separation in 2θ
            # Δ(2θ) ≈ 2(Δλ/λ)tan(θ) where Δλ = λ_Kα2 - λ_Kα1
            delta_lambda = self.ka2 - self.ka1
            avg_lambda = (self.ka1 + self.ka2) / 2
            separation_2theta = 2 * (delta_lambda / avg_lambda) * np.tan(theta_rad)
            separation_deg = np.rad2deg(separation_2theta)

            # Estimate Kα1 component width
            # Simple approximation: assume Gaussian profiles
            # β_observed² ≈ β_Kα1² + (I_Kα2/I_Kα1) * overlap_factor

            # Overlap factor depends on separation relative to peak width
            overlap_ratio = separation_deg / beta_corrected

            if overlap_ratio < 0.5:  # Heavily overlapping
                # Use empirical correction factor
                deconv_factor = 0.85  # Kα1 component is ~85% of observed width
                warnings = ['Severe Kα doublet overlap - deconvolution uncertain']
            elif overlap_ratio < 1.0:  # Moderately overlapping
                deconv_factor = 0.90
                warnings = ['Moderate Kα doublet overlap']
            else:  # Well separated
                deconv_factor = 0.95
                warnings = []

            beta_deconv = beta_corrected * deconv_factor

            doublet_info = {
                'method': 'doublet_deconvolution',
                'correction_applied': True,
                'separation_deg': separation_deg,
                'overlap_ratio': overlap_ratio,
                'deconv_factor': deconv_factor,
                'ka1_wavelength': self.ka1,
                'ka2_wavelength': self.ka2,
                'warnings': warnings
            }

            return beta_deconv, doublet_info

        except Exception as e:
            self.logger.warning(f"Doublet deconvolution failed for peak {peak_number}: {e}")
            return beta_corrected, {
                'method': 'doublet_deconvolution',
                'correction_applied': False,
                'error': str(e),
                'warnings': ['Doublet deconvolution failed, using instrumental-corrected FWHM']
            }

    def _estimate_size_uncertainty(self, size_nm: float, beta_monshi: float,
                                   beta_corrected: float, peak: Dict, doublet_info: Dict) -> float:
        """Estimate size uncertainty including doublet correction uncertainty"""
        try:
            # Base uncertainty from peak fitting
            base_uncertainty = super()._estimate_size_uncertainty(
                size_nm, beta_corrected, peak.get('FWHM_obs', beta_corrected), peak)

            # Additional uncertainty from doublet correction
            if doublet_info.get('correction_applied', False):
                correction_magnitude = doublet_info.get('correction_magnitude', 0)

                # Estimate alpha parameter uncertainty (typically ~10% for well-characterized sources)
                alpha_uncertainty_rel = 0.1

                # Propagate alpha uncertainty to size uncertainty
                doublet_uncertainty_rel = correction_magnitude * alpha_uncertainty_rel
                doublet_uncertainty = size_nm * doublet_uncertainty_rel

                # Combine uncertainties in quadrature
                total_uncertainty = np.sqrt(base_uncertainty ** 2 + doublet_uncertainty ** 2)

                return total_uncertainty
            else:
                return base_uncertainty

        except Exception as e:
            self.logger.debug(f"Monshi size uncertainty estimation failed: {e}")
            return size_nm * 0.15  # Default 15% uncertainty for Monshi method

    def _calculate_quality_metrics(self, results: List[Dict]) -> Dict:
        """Calculate quality metrics specific to Monshi-Scherrer analysis"""
        if not results:
            return {}

        base_metrics = super()._calculate_quality_metrics(results)

        # Monshi-specific metrics
        doublet_corrections = [r['doublet_correction_factor'] for r in results]
        correction_magnitudes = []

        for r in results:
            doublet_info = r.get('doublet_info', {})
            if 'correction_magnitude' in doublet_info:
                correction_magnitudes.append(doublet_info['correction_magnitude'])

        monshi_metrics = {
            'average_doublet_correction': np.mean(doublet_corrections),
            'min_doublet_correction': min(doublet_corrections),
            'large_corrections': sum(1 for cf in doublet_corrections if cf < 0.8),
            'alpha_parameter_used': self.alpha
        }

        if correction_magnitudes:
            monshi_metrics.update({
                'average_correction_magnitude': np.mean(correction_magnitudes),
                'max_correction_magnitude': max(correction_magnitudes)
            })

        # Combine base and Monshi metrics
        base_metrics.update(monshi_metrics)
        return base_metrics

    def _assess_analysis_quality(self, analysis_result: Dict) -> Dict:
        """Enhanced quality assessment for Monshi-Scherrer method"""
        # Start with base assessment
        assessment = {
            'overall_quality': 'unknown',
            'confidence_level': 'unknown',
            'recommendations': [],
            'warnings': []
        }

        try:
            quality_metrics = analysis_result.get('quality_metrics', {})
            doublet_warnings = analysis_result.get('doublet_warnings', [])

            # Check doublet correction impact
            avg_correction = quality_metrics.get('average_doublet_correction', 1.0)
            large_corrections = quality_metrics.get('large_corrections', 0)
            num_peaks = analysis_result.get('num_peaks', 0)

            if avg_correction < 0.7:
                assessment['warnings'].append("Large average doublet corrections - verify alpha parameter")

            if large_corrections > num_peaks * 0.3:
                assessment['warnings'].append("Many peaks require large doublet corrections")

            if doublet_warnings:
                assessment['warnings'].append(f"{len(doublet_warnings)} doublet-related warnings")

            # Wavelength-specific recommendations
            assessment['recommendations'].append(f"Analysis uses Kα1 = {self.ka1:.5f} nm")

            if self.alpha > 3e-3:
                assessment['recommendations'].append("Large alpha parameter - suitable for Cu radiation")
            elif self.alpha < 1.5e-3:
                assessment['recommendations'].append("Small alpha parameter - check source characteristics")

            # Overall quality based on base metrics and Monshi-specific factors
            avg_r2 = quality_metrics.get('average_r_squared', 0)

            if avg_r2 > 0.85 and avg_correction > 0.8 and large_corrections == 0:
                assessment['overall_quality'] = 'excellent'
                assessment['confidence_level'] = 'high'
            elif avg_r2 > 0.75 and avg_correction > 0.7:
                assessment['overall_quality'] = 'good'
                assessment['confidence_level'] = 'moderate_to_high'
            elif avg_r2 > 0.6:
                assessment['overall_quality'] = 'moderate'
                assessment['confidence_level'] = 'moderate'
            else:
                assessment['overall_quality'] = 'poor'
                assessment['confidence_level'] = 'low'

        except Exception as e:
            self.logger.warning(f"Monshi analysis quality assessment failed: {e}")
            assessment['warnings'].append("Quality assessment incomplete")

        return assessment