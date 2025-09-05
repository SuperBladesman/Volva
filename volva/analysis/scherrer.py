"""
Volva XRD Analysis Tool - Scherrer Analysis
Enhanced Scherrer analysis with proper scientific validation
"""

import numpy as np
from typing import List, Dict, Union
import logging

from ..core import AnalysisError, Constants, VolvaConfig
from .base import AnalysisMethod

logger = logging.getLogger(__name__)


class ScherrerAnalyzer(AnalysisMethod):
    """Enhanced Scherrer analysis implementation with scientific robustness"""

    def analyze(self, peaks_data: List[Dict], k_constant: Union[float, str] = None) -> Dict:
        """Enhanced Scherrer analysis with comprehensive validation"""
        self._validate_peaks(peaks_data)

        # Handle K constant with validation
        k_value, k_description = self._process_k_constant(k_constant)

        results = []
        failed_peaks = []
        validation_warnings = []

        for i, peak in enumerate(peaks_data):
            try:
                peak_result = self._analyze_single_peak(peak, k_value, i + 1)
                if peak_result['valid']:
                    results.append(peak_result['data'])
                    validation_warnings.extend(peak_result['warnings'])
                else:
                    failed_peaks.append({
                        'peak_number': i + 1,
                        'reason': peak_result['reason'],
                        '2theta': peak['2theta']
                    })

            except Exception as e:
                self.logger.warning(f"Scherrer calculation failed for peak {i + 1}: {e}")
                failed_peaks.append({
                    'peak_number': i + 1,
                    'reason': f'Calculation error: {str(e)}',
                    '2theta': peak['2theta']
                })
                continue

        if not results:
            raise AnalysisError('No valid peaks for Scherrer analysis')

        # Calculate statistical measures
        sizes = [r['size_nm'] for r in results]
        quality_weights = [r['r_squared'] for r in results]

        # Quality-weighted statistics
        weighted_avg = np.average(sizes, weights=quality_weights)

        analysis_result = {
            'success': True,
            'method': 'Scherrer',
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
            'wavelength': self.wavelength,
            'validation_warnings': validation_warnings,
            'quality_metrics': self._calculate_quality_metrics(results)
        }

        # Overall assessment
        analysis_result['assessment'] = self._assess_analysis_quality(analysis_result)

        return analysis_result

    def _process_k_constant(self, k_constant: Union[float, str, None]) -> tuple:
        """Process and validate K constant"""
        if k_constant is None:
            k_constant = self.config.default_scherrer_k

        if isinstance(k_constant, str):
            if k_constant not in Constants.SCHERRER_CONSTANTS:
                available = list(Constants.SCHERRER_CONSTANTS.keys())
                raise ValueError(f"Unknown Scherrer constant: '{k_constant}'. "
                                 f"Available: {available}")
            k_data = Constants.SCHERRER_CONSTANTS[k_constant]
            k_value = k_data["value"]
            k_description = f"{k_constant} ({k_data['description']})"
        else:
            k_value = float(k_constant)
            if not (0.5 <= k_value <= 2.0):
                self.logger.warning(f"Unusual K constant value: {k_value}")
            k_description = f"Custom value: {k_value}"

        return k_value, k_description

    def _analyze_single_peak(self, peak: Dict, k_value: float, peak_number: int) -> Dict:
        """Analyze a single peak with comprehensive validation"""
        two_theta = peak['2theta']
        fwhm_obs = peak['FWHM_obs']
        fwhm_inst = peak.get('FWHM_inst', self.config.default_instrumental_broadening)

        warnings = []

        # Check Scherrer equation applicability
        scherrer_check = self._check_scherrer_limits(two_theta, fwhm_obs)
        warnings.extend(scherrer_check['warnings'])

        # Correct for instrumental broadening
        beta_corrected = self._get_effective_fwhm(peak)

        # Calculate correction factor for logging
        correction_factor = beta_corrected / fwhm_obs
        if correction_factor < 0.7:
            warnings.append(f"Large instrumental correction ({correction_factor:.1%} retained)")

        # Convert to radians for calculation
        beta_rad = np.deg2rad(beta_corrected)
        theta_rad = np.deg2rad(two_theta / 2.0)

        # Scherrer equation: D = K*λ/(β*cos(θ))
        try:
            cos_theta = np.cos(theta_rad)
            if cos_theta <= 0:
                return {
                    'valid': False,
                    'reason': f'Invalid angle: cos(θ) = {cos_theta:.6f}',
                    'warnings': warnings
                }

            crystallite_size = k_value * self.wavelength / (beta_rad * cos_theta)
            crystallite_size_nm = crystallite_size * 1e9

        except (ZeroDivisionError, ValueError) as e:
            return {
                'valid': False,
                'reason': f'Calculation error: {str(e)}',
                'warnings': warnings
            }

        # Validate result
        if not self._validate_crystallite_size(crystallite_size_nm, f"Peak {peak_number}"):
            return {
                'valid': False,
                'reason': f'Unrealistic size: {crystallite_size_nm:.1f} nm',
                'warnings': warnings
            }

        # Calculate uncertainty if possible
        size_uncertainty = self._estimate_size_uncertainty(
            crystallite_size_nm, beta_corrected, fwhm_obs, peak)

        # Assess individual peak quality
        peak_quality = self._assess_peak_quality(peak, scherrer_check, correction_factor)

        peak_result = {
            'peak_number': peak_number,
            'two_theta': two_theta,
            'size_nm': crystallite_size_nm,
            'size_uncertainty_nm': size_uncertainty,
            'fwhm_obs': fwhm_obs,
            'fwhm_corrected': beta_corrected,
            'fwhm_inst': fwhm_inst,
            'correction_factor': correction_factor,
            'r_squared': peak.get('r_squared', 0.85),
            'theta_rad': theta_rad,
            'beta_rad': beta_rad,
            'cos_theta': cos_theta,
            'peak_quality': peak_quality
        }

        return {
            'valid': True,
            'data': peak_result,
            'warnings': warnings
        }

    def _estimate_size_uncertainty(self, size_nm: float, beta_corrected: float,
                                   fwhm_obs: float, peak: Dict) -> float:
        """Estimate crystallite size uncertainty using error propagation"""
        try:
            # Get fitting uncertainty if available
            fit_info = peak.get('fit_info', {})

            if isinstance(fit_info, dict) and 'fit_errors' in fit_info:
                # Use FWHM fitting error if available
                fwhm_error = None
                for key, error in fit_info['fit_errors'].items():
                    if error is not None and ('sigma' in key.lower() or 'gamma' in key.lower()):
                        # Convert parameter error to FWHM error (approximate)
                        if 'sigma' in key.lower():
                            fwhm_error = 2.355 * error  # Gaussian FWHM = 2.355*sigma
                        elif 'gamma' in key.lower():
                            fwhm_error = 2.0 * error  # Lorentzian FWHM = 2*gamma
                        break

                if fwhm_error is not None and fwhm_error > 0:
                    # Error propagation: δD/D = δβ/β
                    relative_error = fwhm_error / beta_corrected
                    return size_nm * relative_error

            # Fallback: estimate based on peak quality
            r_squared = peak.get('r_squared', 0.85)
            if r_squared > 0.9:
                # High quality: ~2% uncertainty
                return size_nm * 0.02
            elif r_squared > 0.8:
                # Good quality: ~5% uncertainty
                return size_nm * 0.05
            elif r_squared > 0.6:
                # Moderate quality: ~10% uncertainty
                return size_nm * 0.10
            else:
                # Poor quality: ~20% uncertainty
                return size_nm * 0.20

        except Exception as e:
            self.logger.debug(f"Size uncertainty estimation failed: {e}")
            return size_nm * 0.10  # Default 10% uncertainty

    def _assess_peak_quality(self, peak: Dict, scherrer_check: Dict,
                             correction_factor: float) -> str:
        """Assess individual peak quality for Scherrer analysis"""
        r_squared = peak.get('r_squared', 0.5)
        intensity = peak.get('intensity', 0)

        quality_factors = []

        # Fit quality
        if r_squared > 0.9:
            quality_factors.append('excellent_fit')
        elif r_squared > 0.8:
            quality_factors.append('good_fit')
        elif r_squared > 0.6:
            quality_factors.append('moderate_fit')
        else:
            quality_factors.append('poor_fit')

        # Instrumental correction impact
        if correction_factor > 0.9:
            quality_factors.append('minimal_correction')
        elif correction_factor > 0.8:
            quality_factors.append('small_correction')
        elif correction_factor > 0.6:
            quality_factors.append('moderate_correction')
        else:
            quality_factors.append('large_correction')

        # Peak intensity (relative indicator)
        if intensity > 1000:
            quality_factors.append('high_intensity')
        elif intensity > 100:
            quality_factors.append('medium_intensity')
        else:
            quality_factors.append('low_intensity')

        # Overall assessment
        if 'excellent_fit' in quality_factors and 'minimal_correction' in quality_factors:
            return 'excellent'
        elif 'poor_fit' in quality_factors or 'large_correction' in quality_factors:
            return 'poor'
        elif 'good_fit' in quality_factors and correction_factor > 0.7:
            return 'good'
        else:
            return 'moderate'

    def _calculate_quality_metrics(self, results: List[Dict]) -> Dict:
        """Calculate overall quality metrics for the analysis"""
        if not results:
            return {}

        r_squareds = [r['r_squared'] for r in results]
        correction_factors = [r['correction_factor'] for r in results]
        qualities = [r['peak_quality'] for r in results]

        return {
            'average_r_squared': np.mean(r_squareds),
            'min_r_squared': min(r_squareds),
            'average_correction_factor': np.mean(correction_factors),
            'min_correction_factor': min(correction_factors),
            'excellent_peaks': qualities.count('excellent'),
            'good_peaks': qualities.count('good'),
            'moderate_peaks': qualities.count('moderate'),
            'poor_peaks': qualities.count('poor'),
            'size_consistency': self._assess_size_consistency([r['size_nm'] for r in results])
        }

    def _assess_size_consistency(self, sizes: List[float]) -> Dict:
        """Assess consistency of crystallite sizes across peaks"""
        if len(sizes) < 2:
            return {'assessment': 'insufficient_data'}

        cv = np.std(sizes) / np.mean(sizes)  # Coefficient of variation

        if cv < 0.1:
            return {'assessment': 'excellent', 'cv': cv, 'description': 'Very consistent sizes'}
        elif cv < 0.2:
            return {'assessment': 'good', 'cv': cv, 'description': 'Good size consistency'}
        elif cv < 0.4:
            return {'assessment': 'moderate', 'cv': cv, 'description': 'Moderate size variation'}
        else:
            return {'assessment': 'poor', 'cv': cv, 'description': 'Large size variation - check sample homogeneity'}

    def _assess_analysis_quality(self, analysis_result: Dict) -> Dict:
        """Overall assessment of the Scherrer analysis quality"""
        assessment = {
            'overall_quality': 'unknown',
            'confidence_level': 'unknown',
            'recommendations': [],
            'warnings': []
        }

        try:
            num_peaks = analysis_result['num_peaks']
            num_failed = len(analysis_result['failed_peaks'])
            quality_metrics = analysis_result['quality_metrics']

            # Success rate
            success_rate = num_peaks / (num_peaks + num_failed) if (num_peaks + num_failed) > 0 else 0

            # Peak count assessment
            if num_peaks >= 5:
                assessment['recommendations'].append("Good statistical basis with 5+ peaks")
            elif num_peaks >= 3:
                assessment['recommendations'].append("Adequate peak count for analysis")
            else:
                assessment['warnings'].append("Limited peak count - results may be less reliable")

            # Quality assessment
            avg_r2 = quality_metrics.get('average_r_squared', 0)
            excellent_peaks = quality_metrics.get('excellent_peaks', 0)
            poor_peaks = quality_metrics.get('poor_peaks', 0)

            if avg_r2 > 0.9 and excellent_peaks >= num_peaks * 0.7:
                assessment['overall_quality'] = 'excellent'
                assessment['confidence_level'] = 'high'
            elif avg_r2 > 0.8 and poor_peaks <= num_peaks * 0.2:
                assessment['overall_quality'] = 'good'
                assessment['confidence_level'] = 'moderate_to_high'
            elif avg_r2 > 0.6:
                assessment['overall_quality'] = 'moderate'
                assessment['confidence_level'] = 'moderate'
                assessment['warnings'].append("Some peaks have moderate fit quality")
            else:
                assessment['overall_quality'] = 'poor'
                assessment['confidence_level'] = 'low'
                assessment['warnings'].append("Poor average fit quality - results may be unreliable")

            # Size consistency
            size_consistency = quality_metrics.get('size_consistency', {})
            if size_consistency.get('assessment') == 'poor':
                assessment['warnings'].append("Large size variation between peaks")
                assessment['recommendations'].append("Check sample homogeneity and peak assignments")
            elif size_consistency.get('assessment') == 'excellent':
                assessment['recommendations'].append("Excellent size consistency supports reliability")

            # Instrumental correction impact
            avg_correction = quality_metrics.get('average_correction_factor', 1.0)
            if avg_correction < 0.7:
                assessment['warnings'].append("Large instrumental corrections - verify instrument parameters")

        except Exception as e:
            self.logger.warning(f"Analysis quality assessment failed: {e}")
            assessment['warnings'].append("Could not assess analysis quality")

        return assessment