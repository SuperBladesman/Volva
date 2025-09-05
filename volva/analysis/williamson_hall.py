"""
Volva XRD Analysis Tool - Williamson-Hall Analysis
Enhanced Williamson-Hall analysis with isotropy validation and robust statistics
"""

import numpy as np
from typing import List, Dict, Optional, Tuple
import logging
from scipy import stats

from ..core import AnalysisError, InsufficientDataError, VolvaConfig
from .base import StrainAnalysisMethod

logger = logging.getLogger(__name__)


class WilliamsonHallAnalyzer(StrainAnalysisMethod):
    """Enhanced Williamson-Hall analysis with isotropy validation and robust fitting"""

    def analyze(self, peaks_data: List[Dict], k_constant: float = None,
                use_robust_fitting: bool = True, outlier_threshold: float = 2.5) -> Dict:
        """
        Enhanced Williamson-Hall analysis with comprehensive validation

        Parameters:
        - k_constant: Shape factor (default from config)
        - use_robust_fitting: Use robust regression to handle outliers
        - outlier_threshold: Threshold for outlier detection (standard deviations)
        """
        self._validate_peaks(peaks_data, min_peaks=3)

        # Enhanced validation for strain analysis
        strain_validation = self._validate_strain_analysis(peaks_data)
        if not strain_validation['valid']:
            raise InsufficientDataError(f"Strain analysis validation failed: {strain_validation['warnings']}")

        k_value = k_constant if k_constant is not None else self.config.default_wh_k

        if not (0.5 <= k_value <= 2.0):
            self.logger.warning(f"Unusual K constant for W-H analysis: {k_value}")

        # Calculate W-H coordinates
        wh_data = self._calculate_wh_coordinates(peaks_data)

        if len(wh_data['valid_points']) < 3:
            raise InsufficientDataError(f'Need at least 3 valid W-H data points, got {len(wh_data["valid_points"])}')

        # Perform linear regression with outlier handling
        fit_results = self._perform_wh_regression(wh_data, use_robust_fitting, outlier_threshold)

        # Extract physical parameters
        physical_params = self._extract_physical_parameters(fit_results, k_value)

        # Assess isotropy assumption
        isotropy_assessment = self._assess_isotropy_assumption(
            wh_data['x_values'], wh_data['y_values'], fit_results['r_squared'])

        # Calculate uncertainties
        uncertainties = self._calculate_uncertainties(fit_results, k_value, len(wh_data['valid_points']))

        # Compile comprehensive results
        analysis_result = {
            'success': True,
            'method': 'Williamson-Hall',
            'crystallite_size_nm': physical_params['size_nm'],
            'crystallite_size_uncertainty_nm': uncertainties['size_uncertainty'],
            'microstrain_percent': physical_params['strain_percent'],
            'microstrain_uncertainty_percent': uncertainties['strain_uncertainty'],
            'r_squared': fit_results['r_squared'],
            'adjusted_r_squared': fit_results.get('adjusted_r_squared', fit_results['r_squared']),

            # Linear regression details
            'slope': fit_results['slope'],
            'intercept': fit_results['intercept'],
            'slope_error': fit_results.get('slope_error', np.inf),
            'intercept_error': fit_results.get('intercept_error', np.inf),

            # Data and quality metrics
            'x_data': wh_data['x_values'],
            'y_data': wh_data['y_values'],
            'fitted_line': fit_results['fitted_values'],
            'residuals': fit_results['residuals'],
            'peak_info': wh_data['peak_details'],
            'outliers': fit_results.get('outliers', []),
            'num_peaks': len(wh_data['valid_points']),
            'failed_peaks': wh_data['failed_peaks'],

            # Analysis parameters
            'k_constant': k_value,
            'wavelength': self.wavelength,
            'use_robust_fitting': use_robust_fitting,
            'outlier_threshold': outlier_threshold,

            # Quality assessments
            'strain_validation': strain_validation,
            'isotropy_assessment': isotropy_assessment,
            'statistical_tests': self._perform_statistical_tests(fit_results),
            'quality_metrics': self._calculate_wh_quality_metrics(fit_results, wh_data),

            # Overall assessment
            'assessment': {}  # Will be populated below
        }

        # Generate overall assessment
        analysis_result['assessment'] = self._assess_wh_analysis_quality(analysis_result)

        return analysis_result

    def _calculate_wh_coordinates(self, peaks_data: List[Dict]) -> Dict:
        """Calculate Williamson-Hall coordinates with validation"""
        valid_points = []
        failed_peaks = []
        peak_details = []

        for i, peak in enumerate(peaks_data):
            try:
                two_theta = peak['2theta']
                fwhm_obs = peak['FWHM_obs']
                fwhm_inst = peak.get('FWHM_inst', self.config.default_instrumental_broadening)

                # Correct for instrumental broadening
                beta_corrected = self._get_effective_fwhm(peak)
                beta_rad = np.deg2rad(beta_corrected)
                theta_rad = np.deg2rad(two_theta / 2.0)

                # Williamson-Hall coordinates: x = 4sinθ/λ, y = βcosθ/λ
                x_val = 4 * np.sin(theta_rad) / self.wavelength
                y_val = beta_rad * np.cos(theta_rad) / self.wavelength

                # Validate coordinates
                if not (np.isfinite(x_val) and np.isfinite(y_val) and x_val > 0 and y_val > 0):
                    raise ValueError(f"Invalid W-H coordinates: x={x_val}, y={y_val}")

                # Check physical reasonableness
                if x_val > 1e9:  # Unreasonably large scattering vector
                    raise ValueError(f"Unreasonable scattering vector: {x_val:.2e} nm⁻¹")

                valid_points.append(i)
                peak_details.append({
                    'peak_number': i + 1,
                    'two_theta': two_theta,
                    'x_value': x_val,
                    'y_value': y_val,
                    'fwhm_obs': fwhm_obs,
                    'fwhm_corrected': beta_corrected,
                    'fwhm_inst': fwhm_inst,
                    'correction_factor': beta_corrected / fwhm_obs,
                    'theta_rad': theta_rad,
                    'r_squared': peak.get('r_squared', 0.85)
                })

            except Exception as e:
                self.logger.warning(f"W-H coordinate calculation failed for peak {i + 1}: {e}")
                failed_peaks.append({
                    'peak_number': i + 1,
                    'reason': str(e),
                    '2theta': peak.get('2theta', 'unknown')
                })

        if len(valid_points) < 3:
            raise InsufficientDataError(f"Insufficient valid peaks for W-H analysis: {len(valid_points)}")

        x_values = np.array([pd['x_value'] for pd in peak_details])
        y_values = np.array([pd['y_value'] for pd in peak_details])

        return {
            'valid_points': valid_points,
            'failed_peaks': failed_peaks,
            'peak_details': peak_details,
            'x_values': x_values,
            'y_values': y_values
        }

    def _perform_wh_regression(self, wh_data: Dict, use_robust: bool,
                               outlier_threshold: float) -> Dict:
        """Perform Williamson-Hall linear regression with robust statistics"""
        x_data = wh_data['x_values']
        y_data = wh_data['y_values']

        if np.std(x_data) == 0:
            raise AnalysisError("No variation in x-data for linear regression")

        # Standard linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)
        r_squared = r_value ** 2

        # Calculate fitted values and residuals
        y_fitted = slope * x_data + intercept
        residuals = y_data - y_fitted

        # Calculate standard errors
        n = len(x_data)
        x_mean = np.mean(x_data)
        sxx = np.sum((x_data - x_mean) ** 2)

        # Standard error of the estimate
        mse = np.sum(residuals ** 2) / (n - 2)
        s_yx = np.sqrt(mse)

        # Parameter standard errors
        slope_error = s_yx / np.sqrt(sxx) if sxx > 0 else np.inf
        intercept_error = s_yx * np.sqrt(1 / n + x_mean ** 2 / sxx) if sxx > 0 else np.inf

        # Calculate adjusted R²
        adjusted_r_squared = 1 - (1 - r_squared) * (n - 1) / (n - 2) if n > 2 else r_squared

        fit_results = {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_squared,
            'adjusted_r_squared': adjusted_r_squared,
            'p_value': p_value,
            'slope_error': slope_error,
            'intercept_error': intercept_error,
            'fitted_values': y_fitted,
            'residuals': residuals,
            'mse': mse,
            'std_error': s_yx,
            'n_points': n,
            'outliers': []
        }

        # Outlier detection and robust fitting
        if use_robust and n >= 4:
            outlier_results = self._detect_and_handle_outliers(
                x_data, y_data, fit_results, outlier_threshold)
            fit_results.update(outlier_results)

        return fit_results

    def _detect_and_handle_outliers(self, x_data: np.ndarray, y_data: np.ndarray,
                                    initial_fit: Dict, threshold: float) -> Dict:
        """Detect outliers and optionally refit without them"""
        residuals = initial_fit['residuals']
        residual_std = np.std(residuals)

        # Identify outliers using standardized residuals
        outlier_mask = np.abs(residuals) > threshold * residual_std
        outlier_indices = np.where(outlier_mask)[0]

        outlier_info = []
        for idx in outlier_indices:
            outlier_info.append({
                'index': int(idx),
                'x_value': x_data[idx],
                'y_value': y_data[idx],
                'residual': residuals[idx],
                'standardized_residual': residuals[idx] / residual_std
            })

        # If outliers found and enough points remain, perform robust fit
        if len(outlier_indices) > 0 and (len(x_data) - len(outlier_indices)) >= 3:
            self.logger.info(f"Found {len(outlier_indices)} outlier(s), performing robust fit")

            # Create mask for non-outliers
            inlier_mask = ~outlier_mask
            x_robust = x_data[inlier_mask]
            y_robust = y_data[inlier_mask]

            # Refit without outliers
            slope_r, intercept_r, r_value_r, p_value_r, std_err_r = stats.linregress(x_robust, y_robust)

            # Calculate robust statistics
            y_fitted_robust = slope_r * x_data + intercept_r  # Predict for all points
            residuals_robust = y_data - y_fitted_robust

            # Robust fit quality metrics
            n_robust = len(x_robust)
            r_squared_robust = r_value_r ** 2

            # Calculate robust parameter errors
            x_mean_r = np.mean(x_robust)
            sxx_r = np.sum((x_robust - x_mean_r) ** 2)
            mse_r = np.sum((y_robust - slope_r * x_robust - intercept_r) ** 2) / (n_robust - 2)
            s_yx_r = np.sqrt(mse_r)

            slope_error_r = s_yx_r / np.sqrt(sxx_r) if sxx_r > 0 else np.inf
            intercept_error_r = s_yx_r * np.sqrt(1 / n_robust + x_mean_r ** 2 / sxx_r) if sxx_r > 0 else np.inf

            return {
                'robust_fit_used': True,
                'outliers': outlier_info,
                'slope': slope_r,
                'intercept': intercept_r,
                'r_squared': r_squared_robust,
                'slope_error': slope_error_r,
                'intercept_error': intercept_error_r,
                'fitted_values': y_fitted_robust,
                'residuals': residuals_robust,
                'n_points_used': n_robust,
                'n_outliers': len(outlier_indices)
            }
        else:
            return {
                'robust_fit_used': False,
                'outliers': outlier_info
            }

    def _extract_physical_parameters(self, fit_results: Dict, k_value: float) -> Dict:
        """Extract crystallite size and microstrain from W-H fit"""
        slope = fit_results['slope']
        intercept = fit_results['intercept']

        # Williamson-Hall equations:
        # Slope = 4ε (microstrain in decimal)
        # Intercept = Kλ/D (where D is crystallite size)

        # Calculate microstrain
        microstrain_decimal = slope / 4
        microstrain_percent = microstrain_decimal * 100

        # Handle negative strain (possible with noise)
        if microstrain_percent < 0:
            self.logger.warning(f"Negative microstrain calculated ({microstrain_percent:.4f}%), "
                                f"setting to zero - may indicate size-only broadening")
            microstrain_percent = 0.0

        # Calculate crystallite size
        if intercept <= 0:
            self.logger.warning(f"Non-positive intercept ({intercept:.6e}), size calculation unreliable")
            size_nm = float('inf')  # Indicates problematic fit
        else:
            crystallite_size_m = k_value * self.wavelength / intercept
            size_nm = crystallite_size_m * 1e9

        # Validate physical reasonableness
        if size_nm > 0 and not self._validate_crystallite_size(size_nm, "Williamson-Hall"):
            self.logger.warning("Crystallite size outside reasonable range")

        if microstrain_percent > 5.0:  # >5% is very high
            self.logger.warning(f"Very high microstrain ({microstrain_percent:.3f}%) - check data quality")

        return {
            'size_nm': size_nm,
            'strain_percent': microstrain_percent,
            'strain_decimal': microstrain_decimal
        }

    def _calculate_uncertainties(self, fit_results: Dict, k_value: float, n_points: int) -> Dict:
        """Calculate uncertainties using error propagation"""
        slope_error = fit_results.get('slope_error', np.inf)
        intercept_error = fit_results.get('intercept_error', np.inf)
        slope = fit_results['slope']
        intercept = fit_results['intercept']

        # Microstrain uncertainty: δε = δslope / 4
        strain_uncertainty, strain_status = self._calculate_strain_uncertainty(
            slope, slope_error, 'williamson_hall')

        # Crystallite size uncertainty: δD/D = δintercept/intercept
        if intercept > 0 and np.isfinite(intercept_error):
            size_uncertainty_rel = intercept_error / intercept
            size_nm = k_value * self.wavelength * 1e9 / intercept
            size_uncertainty = size_nm * size_uncertainty_rel
        else:
            size_uncertainty = np.inf

        return {
            'strain_uncertainty': strain_uncertainty,
            'strain_uncertainty_status': strain_status,
            'size_uncertainty': size_uncertainty,
            'parameter_correlation': self._calculate_parameter_correlation(fit_results)
        }

    def _calculate_parameter_correlation(self, fit_results: Dict) -> Dict:
        """Calculate correlation between size and strain parameters"""
        try:
            x_data = np.array([fit_results['fitted_values'][i] / fit_results['slope'] - fit_results['intercept'] /
                               fit_results['slope']
                               for i in range(len(fit_results['fitted_values']))])
            x_mean = np.mean(x_data)
            sxx = np.sum((x_data - x_mean) ** 2)
            n = fit_results['n_points']

            # Correlation coefficient between slope and intercept
            correlation_coeff = -x_mean / np.sqrt(sxx / n + x_mean ** 2)

            return {
                'slope_intercept_correlation': correlation_coeff,
                'interpretation': 'negative' if correlation_coeff < -0.5 else 'weak' if abs(
                    correlation_coeff) < 0.5 else 'positive'
            }
        except:
            return {'slope_intercept_correlation': 0.0, 'interpretation': 'unknown'}

    def _perform_statistical_tests(self, fit_results: Dict) -> Dict:
        """Perform additional statistical tests on the fit"""
        residuals = fit_results['residuals']
        n = len(residuals)

        tests = {}

        # Normality test (Shapiro-Wilk if n < 50, otherwise Anderson-Darling)
        if n >= 3:
            try:
                if n <= 50:
                    shapiro_stat, shapiro_p = stats.shapiro(residuals)
                    tests['normality'] = {
                        'test': 'shapiro_wilk',
                        'statistic': shapiro_stat,
                        'p_value': shapiro_p,
                        'normal': shapiro_p > 0.05
                    }
                else:
                    anderson_result = stats.anderson(residuals, dist='norm')
                    tests['normality'] = {
                        'test': 'anderson_darling',
                        'statistic': anderson_result.statistic,
                        'critical_values': anderson_result.critical_values,
                        'normal': anderson_result.statistic < anderson_result.critical_values[2]  # 5% level
                    }
            except:
                tests['normality'] = {'test': 'failed', 'normal': None}

        # Durbin-Watson test for autocorrelation
        if n >= 4:
            try:
                dw_stat = np.sum(np.diff(residuals) ** 2) / np.sum(residuals ** 2)
                tests['autocorrelation'] = {
                    'durbin_watson_statistic': dw_stat,
                    'interpretation': 'positive' if dw_stat < 1.5 else 'none' if dw_stat < 2.5 else 'negative'
                }
            except:
                tests['autocorrelation'] = {'durbin_watson_statistic': None}

        # Homoscedasticity (constant variance) - Breusch-Pagan test approximation
        if n >= 5:
            try:
                fitted_values = fit_results['fitted_values']
                # Simple test: correlation between |residuals| and fitted values
                abs_residuals = np.abs(residuals)
                correlation, p_val = stats.pearsonr(fitted_values, abs_residuals)
                tests['homoscedasticity'] = {
                    'correlation_abs_residuals_fitted': correlation,
                    'p_value': p_val,
                    'homoscedastic': p_val > 0.05
                }
            except:
                tests['homoscedasticity'] = {'homoscedastic': None}

        return tests

    def _calculate_wh_quality_metrics(self, fit_results: Dict, wh_data: Dict) -> Dict:
        """Calculate W-H specific quality metrics"""
        metrics = {
            'fit_quality': 'excellent' if fit_results['r_squared'] > 0.95 else
            'good' if fit_results['r_squared'] > 0.9 else
            'moderate' if fit_results['r_squared'] > 0.8 else 'poor',
            'data_range_coverage': self._assess_data_range_coverage(wh_data['x_values']),
            'peak_quality_distribution': self._assess_peak_quality_distribution(wh_data['peak_details']),
            'regression_diagnostics': self._assess_regression_diagnostics(fit_results)
        }

        return metrics

    def _assess_data_range_coverage(self, x_values: np.ndarray) -> Dict:
        """Assess the coverage of the scattering vector range"""
        x_range = np.max(x_values) - np.min(x_values)
        x_mean = np.mean(x_values)

        return {
            'range_nm_inv': x_range,
            'relative_range': x_range / x_mean if x_mean > 0 else 0,
            'assessment': 'excellent' if x_range > x_mean * 0.8 else
            'good' if x_range > x_mean * 0.5 else
            'moderate' if x_range > x_mean * 0.3 else 'limited'
        }

    def _assess_peak_quality_distribution(self, peak_details: List[Dict]) -> Dict:
        """Assess the quality distribution of peaks used"""
        qualities = [pd['r_squared'] for pd in peak_details]
        correction_factors = [pd['correction_factor'] for pd in peak_details]

        return {
            'average_peak_r_squared': np.mean(qualities),
            'min_peak_r_squared': min(qualities),
            'high_quality_peaks': sum(1 for q in qualities if q > 0.9),
            'average_correction_factor': np.mean(correction_factors),
            'large_corrections': sum(1 for cf in correction_factors if cf < 0.7)
        }

    def _assess_regression_diagnostics(self, fit_results: Dict) -> Dict:
        """Assess regression diagnostic metrics"""
        diagnostics = {
            'leverage_points': 0,  # Points with high influence
            'cook_distance_max': 0.0,  # Maximum Cook's distance
            'condition_number': 1.0  # Multicollinearity measure (always 1 for simple regression)
        }

        # For simple linear regression, we can calculate some basic diagnostics
        n = fit_results['n_points']
        if n > 2:
            # Residual standard error relative to data range
            residuals = fit_results['residuals']
            fitted_values = fit_results['fitted_values']

            diagnostics['residual_se_relative'] = fit_results['std_error'] / np.std(fitted_values)
            diagnostics['max_studentized_residual'] = np.max(np.abs(residuals)) / fit_results['std_error']

        return diagnostics

    def _assess_wh_analysis_quality(self, analysis_result: Dict) -> Dict:
        """Comprehensive quality assessment for W-H analysis"""
        assessment = {
            'overall_quality': 'unknown',
            'confidence_level': 'unknown',
            'scientific_validity': 'unknown',
            'recommendations': [],
            'warnings': [],
            'critical_issues': []
        }

        try:
            r_squared = analysis_result['r_squared']
            isotropy = analysis_result['isotropy_assessment']
            n_peaks = analysis_result['num_peaks']
            outliers = len(analysis_result.get('outliers', []))

            # Fit quality assessment
            if r_squared > 0.95:
                assessment['overall_quality'] = 'excellent'
                assessment['confidence_level'] = 'high'
            elif r_squared > 0.9:
                assessment['overall_quality'] = 'good'
                assessment['confidence_level'] = 'moderate_to_high'
            elif r_squared > 0.8:
                assessment['overall_quality'] = 'moderate'
                assessment['confidence_level'] = 'moderate'
                assessment['warnings'].append(f"Moderate fit quality (R² = {r_squared:.3f})")
            else:
                assessment['overall_quality'] = 'poor'
                assessment['confidence_level'] = 'low'
                assessment['warnings'].append(f"Poor fit quality (R² = {r_squared:.3f})")

            # Isotropy assessment
            if isotropy.get('likely_isotropic', False):
                assessment['scientific_validity'] = 'good'
                if isotropy.get('confidence') == 'very_high':
                    assessment['recommendations'].append("Excellent isotropic fit supports W-H model validity")
            else:
                assessment['scientific_validity'] = 'questionable'
                assessment['critical_issues'].append("Isotropy assumption may be violated")
                assessment['recommendations'].extend(isotropy.get('recommendations', []))

            # Data sufficiency
            if n_peaks >= 6:
                assessment['recommendations'].append("Good statistical basis with 6+ peaks")
            elif n_peaks >= 4:
                assessment['recommendations'].append("Adequate peak count for W-H analysis")
            else:
                assessment['warnings'].append("Limited peak count for strain analysis")

            # Outlier impact
            if outliers > 0:
                outlier_fraction = outliers / (n_peaks + outliers)
                if outlier_fraction > 0.3:
                    assessment['warnings'].append(f"High outlier rate ({outlier_fraction:.1%})")
                else:
                    assessment['recommendations'].append("Robust fitting handled outliers effectively")

            # Physical parameter validation
            size_nm = analysis_result['crystallite_size_nm']
            strain_pct = analysis_result['microstrain_percent']

            if np.isfinite(size_nm) and 1 <= size_nm <= 1000:
                assessment['recommendations'].append(f"Reasonable crystallite size ({size_nm:.0f} nm)")
            else:
                assessment['critical_issues'].append("Unrealistic crystallite size calculation")

            if 0 <= strain_pct <= 2.0:
                assessment['recommendations'].append(f"Reasonable microstrain ({strain_pct:.3f}%)")
            elif strain_pct > 2.0:
                assessment['warnings'].append(f"High microstrain ({strain_pct:.3f}%) - verify results")

        except Exception as e:
            assessment['warnings'].append(f"Quality assessment incomplete: {e}")

        return assessment