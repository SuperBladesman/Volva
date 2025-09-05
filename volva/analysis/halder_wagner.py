"""
Volva XRD Analysis Tool - Halder-Wagner Analysis
Enhanced Halder-Wagner analysis for anisotropic strain with robust validation
"""

import numpy as np
from typing import List, Dict, Optional, Tuple
import logging
from scipy import stats
from scipy.optimize import curve_fit

from ..core import AnalysisError, InsufficientDataError, VolvaConfig
from .base import StrainAnalysisMethod

logger = logging.getLogger(__name__)


class HalderWagnerAnalyzer(StrainAnalysisMethod):
    """Enhanced Halder-Wagner analysis for anisotropic and non-uniform strain"""

    def analyze(self, peaks_data: List[Dict], a_parameter: float = None,
                k_constant: float = None, use_robust_fitting: bool = True,
                compare_with_wh: bool = True) -> Dict:
        """
        Enhanced Halder-Wagner analysis with anisotropic strain modeling

        Parameters:
        - a_parameter: Anisotropy parameter (default from config)
        - k_constant: Shape factor (default from config)
        - use_robust_fitting: Use robust regression for outlier handling
        - compare_with_wh: Compare with Williamson-Hall for model selection
        """
        self._validate_peaks(peaks_data, min_peaks=3)

        # Enhanced validation for strain analysis
        strain_validation = self._validate_strain_analysis(peaks_data)
        if not strain_validation['valid']:
            raise InsufficientDataError(f"Strain analysis validation failed: {strain_validation['warnings']}")

        # Set parameters
        a_value = a_parameter if a_parameter is not None else self.config.default_hw_a
        k_value = k_constant if k_constant is not None else self.config.default_hw_k

        # Validate parameters
        if not (1.0 <= a_value <= 4.0):
            self.logger.warning(f"Unusual 'a' parameter for H-W analysis: {a_value}")

        # Calculate H-W coordinates
        hw_data = self._calculate_hw_coordinates(peaks_data)

        if len(hw_data['valid_points']) < 3:
            raise InsufficientDataError(f'Need at least 3 valid H-W data points, got {len(hw_data["valid_points"])}')

        # Perform quadratic regression
        fit_results = self._perform_hw_regression(hw_data, use_robust_fitting)

        # Extract physical parameters
        physical_params = self._extract_physical_parameters(fit_results, k_value, a_value)

        # Compare with Williamson-Hall if requested
        wh_comparison = {}
        if compare_with_wh and len(hw_data['valid_points']) >= 4:
            wh_comparison = self._compare_with_williamson_hall(hw_data)

        # Assess strain anisotropy
        anisotropy_assessment = self._assess_strain_anisotropy(
            hw_data, fit_results, wh_comparison)

        # Calculate uncertainties
        uncertainties = self._calculate_uncertainties(fit_results, k_value, a_value)

        # Compile comprehensive results
        analysis_result = {
            'success': True,
            'method': 'Halder-Wagner',
            'crystallite_size_nm': physical_params['size_nm'],
            'crystallite_size_uncertainty_nm': uncertainties['size_uncertainty'],
            'microstrain_percent': physical_params['strain_percent'],
            'microstrain_uncertainty_percent': uncertainties['strain_uncertainty'],
            'r_squared': fit_results['r_squared'],
            'adjusted_r_squared': fit_results.get('adjusted_r_squared', fit_results['r_squared']),

            # Quadratic regression details
            'slope': fit_results['slope'],
            'intercept': fit_results['intercept'],
            'quadratic_term': fit_results.get('quadratic_coeff', 0.0),
            'slope_error': fit_results.get('slope_error', np.inf),
            'intercept_error': fit_results.get('intercept_error', np.inf),

            # Data and quality metrics
            'x_data': hw_data['x_values'],
            'y_data': hw_data['y_values'],
            'fitted_curve': fit_results['fitted_values'],
            'residuals': fit_results['residuals'],
            'peak_info': hw_data['peak_details'],
            'outliers': fit_results.get('outliers', []),
            'num_peaks': len(hw_data['valid_points']),
            'failed_peaks': hw_data['failed_peaks'],

            # Analysis parameters
            'a_parameter': a_value,
            'k_constant': k_value,
            'wavelength': self.wavelength,
            'use_robust_fitting': use_robust_fitting,

            # Quality assessments
            'strain_validation': strain_validation,
            'anisotropy_assessment': anisotropy_assessment,
            'wh_comparison': wh_comparison,
            'statistical_tests': self._perform_statistical_tests(fit_results),
            'quality_metrics': self._calculate_hw_quality_metrics(fit_results, hw_data),

            # Overall assessment
            'assessment': {}  # Will be populated below
        }

        # Generate overall assessment
        analysis_result['assessment'] = self._assess_hw_analysis_quality(analysis_result)

        return analysis_result

    def _calculate_hw_coordinates(self, peaks_data: List[Dict]) -> Dict:
        """Calculate Halder-Wagner coordinates with validation"""
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

                # Halder-Wagner coordinates: x = (4sinθ/λ)², y = (βcosθ/λ)²
                sin_theta_over_lambda = np.sin(theta_rad) / self.wavelength
                cos_theta = np.cos(theta_rad)

                x_val = (4 * sin_theta_over_lambda) ** 2
                y_val = (beta_rad * cos_theta / self.wavelength) ** 2

                # Validate coordinates
                if not (np.isfinite(x_val) and np.isfinite(y_val) and x_val > 0 and y_val > 0):
                    raise ValueError(f"Invalid H-W coordinates: x={x_val}, y={y_val}")

                # Check physical reasonableness
                if x_val > 1e18:  # Unreasonably large squared scattering vector
                    raise ValueError(f"Unreasonable squared scattering vector: {x_val:.2e} nm⁻²")

                valid_points.append(i)
                peak_details.append({
                    'peak_number': i + 1,
                    'two_theta': two_theta,
                    'x_value': x_val,
                    'y_value': y_val,
                    'x_linear': 4 * sin_theta_over_lambda,  # For W-H comparison
                    'y_linear': beta_rad * cos_theta / self.wavelength,
                    'fwhm_obs': fwhm_obs,
                    'fwhm_corrected': beta_corrected,
                    'fwhm_inst': fwhm_inst,
                    'correction_factor': beta_corrected / fwhm_obs,
                    'theta_rad': theta_rad,
                    'r_squared': peak.get('r_squared', 0.85)
                })

            except Exception as e:
                self.logger.warning(f"H-W coordinate calculation failed for peak {i + 1}: {e}")
                failed_peaks.append({
                    'peak_number': i + 1,
                    'reason': str(e),
                    '2theta': peak.get('2theta', 'unknown')
                })

        if len(valid_points) < 3:
            raise InsufficientDataError(f"Insufficient valid peaks for H-W analysis: {len(valid_points)}")

        x_values = np.array([pd['x_value'] for pd in peak_details])
        y_values = np.array([pd['y_value'] for pd in peak_details])

        return {
            'valid_points': valid_points,
            'failed_peaks': failed_peaks,
            'peak_details': peak_details,
            'x_values': x_values,
            'y_values': y_values
        }

    def _perform_hw_regression(self, hw_data: Dict, use_robust: bool) -> Dict:
        """Perform Halder-Wagner quadratic/linear regression"""
        x_data = hw_data['x_values']
        y_data = hw_data['y_values']
        n = len(x_data)

        if np.std(x_data) == 0:
            raise AnalysisError("No variation in x-data for regression")

        # Try both linear and quadratic fits to determine best model
        linear_fit = self._fit_linear_model(x_data, y_data)

        fit_results = linear_fit.copy()
        fit_results['model_type'] = 'linear'

        # Try quadratic fit if we have enough points
        if n >= 4:
            try:
                quadratic_fit = self._fit_quadratic_model(x_data, y_data)

                # Compare models using F-test and information criteria
                model_comparison = self._compare_linear_quadratic(linear_fit, quadratic_fit, n)

                if model_comparison['preferred'] == 'quadratic':
                    fit_results = quadratic_fit.copy()
                    fit_results['model_type'] = 'quadratic'
                    fit_results['model_comparison'] = model_comparison
                    self.logger.info("Quadratic model preferred over linear for H-W analysis")
                else:
                    fit_results['model_comparison'] = model_comparison
                    self.logger.info("Linear model adequate for H-W analysis")

            except Exception as e:
                self.logger.warning(f"Quadratic fitting failed: {e}, using linear model")

        # Handle outliers if requested
        if use_robust and n >= 4:
            outlier_results = self._detect_and_handle_outliers_hw(
                x_data, y_data, fit_results)
            fit_results.update(outlier_results)

        return fit_results

    def _fit_linear_model(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict:
        """Fit linear model to H-W data"""
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)
        r_squared = r_value ** 2

        # Calculate fitted values and residuals
        y_fitted = slope * x_data + intercept
        residuals = y_data - y_fitted

        # Calculate standard errors
        n = len(x_data)
        x_mean = np.mean(x_data)
        sxx = np.sum((x_data - x_mean) ** 2)

        mse = np.sum(residuals ** 2) / (n - 2)
        s_yx = np.sqrt(mse)

        slope_error = s_yx / np.sqrt(sxx) if sxx > 0 else np.inf
        intercept_error = s_yx * np.sqrt(1 / n + x_mean ** 2 / sxx) if sxx > 0 else np.inf

        return {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_squared,
            'adjusted_r_squared': 1 - (1 - r_squared) * (n - 1) / (n - 2) if n > 2 else r_squared,
            'p_value': p_value,
            'slope_error': slope_error,
            'intercept_error': intercept_error,
            'fitted_values': y_fitted,
            'residuals': residuals,
            'mse': mse,
            'aic': n * np.log(mse) + 2 * 2,  # 2 parameters
            'bic': n * np.log(mse) + np.log(n) * 2,
            'n_points': n
        }

    def _fit_quadratic_model(self, x_data: np.ndarray, y_data: np.ndarray) -> Dict:
        """Fit quadratic model to H-W data"""
        # Quadratic model: y = ax² + bx + c
        coeffs = np.polyfit(x_data, y_data, 2)
        a_coeff, b_coeff, c_coeff = coeffs

        y_fitted = np.polyval(coeffs, x_data)
        residuals = y_data - y_fitted

        # Calculate R²
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        n = len(x_data)
        mse = ss_res / (n - 3)  # 3 parameters

        # For H-W, we interpret: intercept = c_coeff, slope = b_coeff
        return {
            'slope': b_coeff,
            'intercept': c_coeff,
            'quadratic_coeff': a_coeff,
            'r_squared': r_squared,
            'adjusted_r_squared': 1 - (1 - r_squared) * (n - 1) / (n - 3) if n > 3 else r_squared,
            'fitted_values': y_fitted,
            'residuals': residuals,
            'mse': mse,
            'aic': n * np.log(mse) + 2 * 3,  # 3 parameters
            'bic': n * np.log(mse) + np.log(n) * 3,
            'n_points': n,
            'coefficients': coeffs
        }

    def _compare_linear_quadratic(self, linear_fit: Dict, quadratic_fit: Dict, n: int) -> Dict:
        """Compare linear vs quadratic models using statistical tests"""
        # F-test for nested models
        ss_res_linear = np.sum(linear_fit['residuals'] ** 2)
        ss_res_quad = np.sum(quadratic_fit['residuals'] ** 2)

        f_stat = ((ss_res_linear - ss_res_quad) / 1) / (ss_res_quad / (n - 3))
        f_p_value = 1 - stats.f.cdf(f_stat, 1, n - 3)

        # Information criteria comparison
        delta_aic = quadratic_fit['aic'] - linear_fit['aic']
        delta_bic = quadratic_fit['bic'] - linear_fit['bic']

        # Decision logic
        if f_p_value < 0.05 and delta_aic < -2:  # Significant improvement and strong evidence
            preferred = 'quadratic'
            confidence = 'high'
        elif f_p_value < 0.05:  # Significant but not strong evidence
            preferred = 'quadratic'
            confidence = 'moderate'
        else:  # No significant improvement
            preferred = 'linear'
            confidence = 'high' if delta_aic > 2 else 'moderate'

        return {
            'preferred': preferred,
            'confidence': confidence,
            'f_statistic': f_stat,
            'f_p_value': f_p_value,
            'delta_aic': delta_aic,
            'delta_bic': delta_bic,
            'r_squared_improvement': quadratic_fit['r_squared'] - linear_fit['r_squared']
        }

    def _detect_and_handle_outliers_hw(self, x_data: np.ndarray, y_data: np.ndarray,
                                       initial_fit: Dict) -> Dict:
        """Detect and handle outliers for H-W analysis"""
        residuals = initial_fit['residuals']
        residual_std = np.std(residuals)

        # Use Cook's distance-like measure for influence
        n = len(x_data)
        leverage = self._calculate_leverage_hw(x_data, initial_fit.get('model_type', 'linear'))

        # Standardized residuals
        std_residuals = residuals / residual_std

        # Combined outlier score (residual + leverage)
        outlier_scores = np.abs(std_residuals) + 2 * leverage
        threshold = 2.5  # Adjusted threshold

        outlier_mask = outlier_scores > threshold
        outlier_indices = np.where(outlier_mask)[0]

        outlier_info = []
        for idx in outlier_indices:
            outlier_info.append({
                'index': int(idx),
                'x_value': x_data[idx],
                'y_value': y_data[idx],
                'residual': residuals[idx],
                'standardized_residual': std_residuals[idx],
                'leverage': leverage[idx],
                'outlier_score': outlier_scores[idx]
            })

        return {
            'outliers': outlier_info,
            'robust_fit_attempted': len(outlier_indices) > 0
        }

    def _calculate_leverage_hw(self, x_data: np.ndarray, model_type: str) -> np.ndarray:
        """Calculate leverage values for H-W regression"""
        try:
            if model_type == 'quadratic':
                # Design matrix for quadratic: [1, x, x²]
                X = np.column_stack([np.ones(len(x_data)), x_data, x_data ** 2])
            else:
                # Design matrix for linear: [1, x]
                X = np.column_stack([np.ones(len(x_data)), x_data])

            # Leverage = diagonal of hat matrix H = X(X'X)⁻¹X'
            try:
                XtX_inv = np.linalg.inv(X.T @ X)
                leverage = np.sum(X * (X @ XtX_inv), axis=1)
                return leverage
            except np.linalg.LinAlgError:
                # Fallback to simple distance measure
                x_mean = np.mean(x_data)
                x_std = np.std(x_data)
                return np.abs(x_data - x_mean) / x_std if x_std > 0 else np.zeros(len(x_data))

        except Exception as e:
            self.logger.debug(f"Leverage calculation failed: {e}")
            return np.zeros(len(x_data))

    def _extract_physical_parameters(self, fit_results: Dict, k_value: float, a_value: float) -> Dict:
        """Extract crystallite size and microstrain from H-W fit"""
        slope = fit_results['slope']
        intercept = fit_results['intercept']

        # Halder-Wagner equations:
        # Slope = ε²/a (where ε is microstrain in decimal, a is anisotropy parameter)
        # Intercept = (Kλ/D)² (where D is crystallite size)

        # Calculate microstrain
        if slope < 0:
            self.logger.warning(f"Negative slope in H-W analysis ({slope:.6e}), "
                                f"setting strain to zero - indicates size-only broadening")
            microstrain_decimal = 0.0
        else:
            microstrain_decimal = np.sqrt(slope * a_value)

        microstrain_percent = microstrain_decimal * 100

        # Calculate crystallite size
        if intercept <= 0:
            self.logger.warning(f"Non-positive intercept ({intercept:.6e}), size calculation unreliable")
            size_nm = float('inf')
        else:
            # D = Kλ/√intercept
            crystallite_size_m = k_value * self.wavelength / np.sqrt(intercept)
            size_nm = crystallite_size_m * 1e9

        # Validate physical reasonableness
        if size_nm > 0 and not self._validate_crystallite_size(size_nm, "Halder-Wagner"):
            self.logger.warning("H-W crystallite size outside reasonable range")

        if microstrain_percent > 5.0:
            self.logger.warning(f"Very high H-W microstrain ({microstrain_percent:.3f}%)")

        return {
            'size_nm': size_nm,
            'strain_percent': microstrain_percent,
            'strain_decimal': microstrain_decimal
        }

    def _compare_with_williamson_hall(self, hw_data: Dict) -> Dict:
        """Compare H-W results with equivalent W-H analysis"""
        try:
            # Extract linear coordinates from H-W data
            x_linear = np.array([pd['x_linear'] for pd in hw_data['peak_details']])
            y_linear = np.array([pd['y_linear'] for pd in hw_data['peak_details']])

            # Perform W-H linear regression
            wh_fit = self._fit_linear_model(x_linear, y_linear)

            # Extract W-H parameters
            wh_slope = wh_fit['slope']
            wh_intercept = wh_fit['intercept']
            wh_r_squared = wh_fit['r_squared']

            # Calculate W-H physical parameters
            wh_strain_percent = (wh_slope / 4) * 100
            wh_size_nm = (
                                     self.config.default_wh_k * self.wavelength * 1e9) / wh_intercept if wh_intercept > 0 else np.inf

            return {
                'comparison_performed': True,
                'wh_r_squared': wh_r_squared,
                'wh_strain_percent': max(0, wh_strain_percent),  # Ensure non-negative
                'wh_size_nm': wh_size_nm,
                'wh_slope': wh_slope,
                'wh_intercept': wh_intercept,
                'method_agreement': self._assess_method_agreement(wh_strain_percent, wh_size_nm)
            }

        except Exception as e:
            self.logger.warning(f"W-H comparison failed: {e}")
            return {'comparison_performed': False, 'error': str(e)}

    def _assess_method_agreement(self, wh_strain: float, wh_size: float) -> str:
        """Assess agreement between H-W and W-H methods (placeholder for now)"""
        # This will be used when we have H-W results to compare
        # For now, return basic assessment
        if wh_strain < 0:
            return 'poor_wh_negative_strain'
        elif wh_strain > 3.0:
            return 'high_strain_check_anisotropy'
        elif np.isinf(wh_size) or wh_size <= 0:
            return 'poor_wh_size_calculation'
        else:
            return 'reasonable_wh_results'

    def _assess_strain_anisotropy(self, hw_data: Dict, fit_results: Dict,
                                  wh_comparison: Dict) -> Dict:
        """Assess whether strain is anisotropic"""
        assessment = {
            'likely_anisotropic': False,
            'evidence': [],
            'confidence': 'low',
            'recommendations': []
        }

        try:
            model_type = fit_results.get('model_type', 'linear')
            r_squared_hw = fit_results['r_squared']

            # Evidence for anisotropy
            evidence_score = 0

            # 1. Model preference
            if model_type == 'quadratic':
                model_comp = fit_results.get('model_comparison', {})
                if model_comp.get('confidence') == 'high':
                    evidence_score += 3
                    assessment['evidence'].append('Quadratic model significantly better than linear')
                elif model_comp.get('confidence') == 'moderate':
                    evidence_score += 2
                    assessment['evidence'].append('Quadratic model moderately better than linear')

            # 2. H-W vs W-H comparison
            if wh_comparison.get('comparison_performed', False):
                wh_r_squared = wh_comparison.get('wh_r_squared', 0)
                r_squared_improvement = r_squared_hw - wh_r_squared

                if r_squared_improvement > 0.05:  # 5% improvement
                    evidence_score += 2
                    assessment['evidence'].append(f'H-W fit much better than W-H (ΔR² = {r_squared_improvement:.3f})')
                elif r_squared_improvement > 0.02:  # 2% improvement
                    evidence_score += 1
                    assessment['evidence'].append(f'H-W fit better than W-H (ΔR² = {r_squared_improvement:.3f})')

                # Check for W-H problems
                wh_agreement = wh_comparison.get('method_agreement', '')
                if 'negative_strain' in wh_agreement:
                    evidence_score += 2
                    assessment['evidence'].append('W-H gives negative strain, H-W resolves this')

            # 3. Residual analysis
            residuals = fit_results.get('residuals', [])
            if len(residuals) > 3:
                # Check for systematic patterns in residuals
                residual_trend = np.corrcoef(hw_data['x_values'], residuals)[0, 1]
                if abs(residual_trend) > 0.5:  # Strong correlation
                    evidence_score += 1
                    assessment['evidence'].append('Systematic residual patterns suggest non-linear relationship')

            # 4. R² quality
            if r_squared_hw > 0.95:
                evidence_score += 1
                assessment['evidence'].append('Excellent H-W fit quality')

            # Overall assessment
            if evidence_score >= 4:
                assessment['likely_anisotropic'] = True
                assessment['confidence'] = 'high'
                assessment['recommendations'].append('Strong evidence for anisotropic strain - H-W method appropriate')
            elif evidence_score >= 2:
                assessment['likely_anisotropic'] = True
                assessment['confidence'] = 'moderate'
                assessment['recommendations'].append('Moderate evidence for anisotropic strain')
            else:
                assessment['confidence'] = 'low'
                assessment['recommendations'].append(
                    'Limited evidence for anisotropy - consider both W-H and H-W results')

            assessment['evidence_score'] = evidence_score

        except Exception as e:
            self.logger.warning(f"Anisotropy assessment failed: {e}")
            assessment['error'] = str(e)

        return assessment

    def _calculate_uncertainties(self, fit_results: Dict, k_value: float, a_value: float) -> Dict:
        """Calculate uncertainties for H-W parameters"""
        slope_error = fit_results.get('slope_error', np.inf)
        intercept_error = fit_results.get('intercept_error', np.inf)
        slope = fit_results['slope']
        intercept = fit_results['intercept']

        # Microstrain uncertainty: ε = √(slope × a), so δε/ε = (1/2) × δslope/slope
        if slope > 0 and np.isfinite(slope_error):
            strain_uncertainty_rel = 0.5 * slope_error / slope
            strain_decimal = np.sqrt(slope * a_value)
            strain_uncertainty = strain_decimal * strain_uncertainty_rel * 100  # Convert to percentage
        else:
            strain_uncertainty = np.inf

        # Size uncertainty: D = Kλ/√intercept, so δD/D = (1/2) × δintercept/intercept
        if intercept > 0 and np.isfinite(intercept_error):
            size_uncertainty_rel = 0.5 * intercept_error / intercept
            size_nm = k_value * self.wavelength * 1e9 / np.sqrt(intercept)
            size_uncertainty = size_nm * size_uncertainty_rel
        else:
            size_uncertainty = np.inf

        return {
            'strain_uncertainty': strain_uncertainty,
            'size_uncertainty': size_uncertainty,
            'parameter_correlation': self._calculate_parameter_correlation_hw(fit_results)
        }

    def _calculate_parameter_correlation_hw(self, fit_results: Dict) -> Dict:
        """Calculate parameter correlations for H-W method"""
        try:
            model_type = fit_results.get('model_type', 'linear')
            # For simple cases, return basic info
            return {
                'model_type': model_type,
                'correlation_available': False,
                'note': 'Detailed correlation analysis not implemented'
            }
        except:
            return {'correlation_available': False}

    def _perform_statistical_tests(self, fit_results: Dict) -> Dict:
        """Perform statistical tests for H-W fit"""
        residuals = fit_results['residuals']
        n = len(residuals)

        tests = {}

        # Normality test
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
                    tests['normality'] = {'test': 'skipped', 'reason': 'too_many_points'}
            except:
                tests['normality'] = {'test': 'failed'}

        # Randomness test for residuals
        if n >= 4:
            try:
                # Runs test for randomness
                median_residual = np.median(residuals)
                runs, n1, n2 = self._runs_test(residuals > median_residual)
                expected_runs = (2 * n1 * n2) / (n1 + n2) + 1
                variance_runs = (2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / ((n1 + n2) ** 2 * (n1 + n2 - 1))

                if variance_runs > 0:
                    z_score = (runs - expected_runs) / np.sqrt(variance_runs)
                    p_value_runs = 2 * (1 - stats.norm.cdf(abs(z_score)))
                    tests['randomness'] = {
                        'runs': runs,
                        'expected_runs': expected_runs,
                        'z_score': z_score,
                        'p_value': p_value_runs,
                        'random': p_value_runs > 0.05
                    }
            except:
                tests['randomness'] = {'test': 'failed'}

        return tests

    def _runs_test(self, binary_sequence: np.ndarray) -> Tuple[int, int, int]:
        """Perform runs test on binary sequence"""
        n1 = np.sum(binary_sequence)
        n2 = len(binary_sequence) - n1

        # Count runs
        runs = 1
        for i in range(1, len(binary_sequence)):
            if binary_sequence[i] != binary_sequence[i - 1]:
                runs += 1

        return runs, n1, n2

    def _calculate_hw_quality_metrics(self, fit_results: Dict, hw_data: Dict) -> Dict:
        """Calculate H-W specific quality metrics"""
        metrics = {
            'fit_quality': 'excellent' if fit_results['r_squared'] > 0.95 else
            'good' if fit_results['r_squared'] > 0.9 else
            'moderate' if fit_results['r_squared'] > 0.8 else 'poor',
            'model_type': fit_results.get('model_type', 'linear'),
            'data_range_coverage': self._assess_hw_data_range_coverage(hw_data['x_values']),
            'peak_quality_distribution': self._assess_peak_quality_distribution_hw(hw_data['peak_details']),
            'regression_diagnostics': self._assess_hw_regression_diagnostics(fit_results)
        }

        # Model-specific metrics
        if fit_results.get('model_type') == 'quadratic':
            model_comp = fit_results.get('model_comparison', {})
            metrics['quadratic_improvement'] = model_comp.get('r_squared_improvement', 0)
            metrics['quadratic_significance'] = model_comp.get('f_p_value', 1.0)

        return metrics

    def _assess_hw_data_range_coverage(self, x_values: np.ndarray) -> Dict:
        """Assess coverage of squared scattering vector range"""
        x_range = np.max(x_values) - np.min(x_values)
        x_mean = np.mean(x_values)

        return {
            'range_nm_inv_squared': x_range,
            'relative_range': x_range / x_mean if x_mean > 0 else 0,
            'assessment': 'excellent' if x_range > x_mean * 0.8 else
            'good' if x_range > x_mean * 0.5 else
            'moderate' if x_range > x_mean * 0.3 else 'limited'
        }

    def _assess_peak_quality_distribution_hw(self, peak_details: List[Dict]) -> Dict:
        """Assess quality distribution for H-W analysis"""
        qualities = [pd['r_squared'] for pd in peak_details]
        correction_factors = [pd['correction_factor'] for pd in peak_details]

        return {
            'average_peak_r_squared': np.mean(qualities),
            'min_peak_r_squared': min(qualities),
            'high_quality_peaks': sum(1 for q in qualities if q > 0.9),
            'average_correction_factor': np.mean(correction_factors),
            'large_corrections': sum(1 for cf in correction_factors if cf < 0.7),
            'peak_count_adequacy': 'excellent' if len(qualities) >= 6 else
            'good' if len(qualities) >= 4 else
            'minimal' if len(qualities) >= 3 else 'insufficient'
        }

    def _assess_hw_regression_diagnostics(self, fit_results: Dict) -> Dict:
        """Assess H-W regression diagnostics"""
        diagnostics = {
            'model_type': fit_results.get('model_type', 'linear'),
            'r_squared': fit_results['r_squared'],
            'adjusted_r_squared': fit_results.get('adjusted_r_squared', fit_results['r_squared'])
        }

        # Information criteria if available
        if 'aic' in fit_results:
            diagnostics['aic'] = fit_results['aic']
        if 'bic' in fit_results:
            diagnostics['bic'] = fit_results['bic']

        # Model comparison metrics
        if 'model_comparison' in fit_results:
            comp = fit_results['model_comparison']
            diagnostics['model_selection'] = {
                'preferred_model': comp.get('preferred', 'unknown'),
                'selection_confidence': comp.get('confidence', 'unknown'),
                'f_test_p_value': comp.get('f_p_value', 1.0)
            }

        return diagnostics

    def _assess_hw_analysis_quality(self, analysis_result: Dict) -> Dict:
        """Comprehensive quality assessment for H-W analysis"""
        assessment = {
            'overall_quality': 'unknown',
            'confidence_level': 'unknown',
            'scientific_validity': 'unknown',
            'model_appropriateness': 'unknown',
            'recommendations': [],
            'warnings': [],
            'critical_issues': []
        }

        try:
            r_squared = analysis_result['r_squared']
            anisotropy = analysis_result['anisotropy_assessment']
            n_peaks = analysis_result['num_peaks']
            model_type = analysis_result.get('quality_metrics', {}).get('model_type', 'linear')

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

            # Model appropriateness
            if anisotropy.get('likely_anisotropic', False):
                if anisotropy.get('confidence') in ['high', 'moderate']:
                    assessment['model_appropriateness'] = 'good'
                    assessment['scientific_validity'] = 'good'
                    assessment['recommendations'].append("H-W model appropriate for anisotropic strain")

                    if model_type == 'quadratic':
                        assessment['recommendations'].append(
                            "Quadratic model successfully captures anisotropic behavior")
                else:
                    assessment['model_appropriateness'] = 'uncertain'
                    assessment['warnings'].append("Uncertain evidence for anisotropy")
            else:
                assessment['model_appropriateness'] = 'questionable'
                assessment['warnings'].append("Limited evidence for anisotropic strain - W-H may be sufficient")
                assessment['recommendations'].append("Compare with Williamson-Hall results")

            # W-H comparison assessment
            wh_comparison = analysis_result.get('wh_comparison', {})
            if wh_comparison.get('comparison_performed', False):
                wh_agreement = wh_comparison.get('method_agreement', '')
                if 'negative_strain' in wh_agreement:
                    assessment['recommendations'].append("H-W resolves negative strain issue from W-H analysis")
                    assessment['scientific_validity'] = 'good'
                elif 'reasonable' in wh_agreement:
                    assessment['recommendations'].append("Both H-W and W-H give reasonable results")

            # Data sufficiency
            if n_peaks >= 6:
                assessment['recommendations'].append("Excellent statistical basis with 6+ peaks")
            elif n_peaks >= 4:
                assessment['recommendations'].append("Good peak count for H-W analysis")
            else:
                assessment['warnings'].append("Limited peak count for strain analysis")

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
                assessment['warnings'].append(f"High microstrain ({strain_pct:.3f}%) - verify anisotropic model")

            # Overall scientific validity
            if (assessment['model_appropriateness'] == 'good' and
                    assessment['overall_quality'] in ['excellent', 'good']):
                assessment['scientific_validity'] = 'good'
            elif assessment['model_appropriateness'] == 'uncertain':
                assessment['scientific_validity'] = 'uncertain'
            else:
                assessment['scientific_validity'] = 'questionable'

        except Exception as e:
            assessment['warnings'].append(f"Quality assessment incomplete: {e}")

        return assessment