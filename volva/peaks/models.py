"""
Volva XRD Analysis Tool - Peak Models
Mathematical models for XRD peak fitting
"""

import numpy as np
import logging
from typing import Tuple, Optional, Dict, Any
from ..core import FittingError, LMFIT_AVAILABLE

logger = logging.getLogger(__name__)

if LMFIT_AVAILABLE:
    from lmfit import models, Model

from typing import NamedTuple

class PeakFitResult(NamedTuple):
    """Result of peak fitting operation."""
    fwhm: float
    r_squared: float
    fit_info: Dict[str, Any]
    center: Optional[float] = None
    amplitude: Optional[float] = None

class PeakModel:
    """Peak models for fitting XRD peaks"""

    @staticmethod
    def gaussian(x: np.ndarray, amplitude: float, center: float, sigma: float) -> np.ndarray:
        """Gaussian peak model"""
        return amplitude * np.exp(-(x - center) ** 2 / (2 * sigma ** 2))

    @staticmethod
    def lorentzian(x: np.ndarray, amplitude: float, center: float, gamma: float) -> np.ndarray:
        """Lorentzian peak model"""
        return amplitude * gamma ** 2 / ((x - center) ** 2 + gamma ** 2)

    @staticmethod
    def voigt_simple(x: np.ndarray, amplitude: float, center: float,
                     sigma: float, gamma: float) -> np.ndarray:
        """Simplified Voigt profile (Gaussian + Lorentzian convolution approximation)"""
        # Simple approximation of Voigt profile
        gauss = np.exp(-(x - center) ** 2 / (2 * sigma ** 2))
        lorentz = gamma ** 2 / ((x - center) ** 2 + gamma ** 2)
        return amplitude * (0.5 * gauss + 0.5 * lorentz)

    @staticmethod
    def pseudo_voigt(x: np.ndarray, amplitude: float, center: float,
                     sigma: float, fraction: float) -> np.ndarray:
        """Pseudo-Voigt profile (linear combination of Gaussian and Lorentzian)"""
        gaussian_part = np.exp(-(x - center) ** 2 / (2 * sigma ** 2))
        lorentzian_part = 1 / (1 + ((x - center) / sigma) ** 2)
        return amplitude * (fraction * lorentzian_part + (1 - fraction) * gaussian_part)


class PeakFitter:
    """Advanced peak fitting with multiple models"""

    def __init__(self, use_lmfit: bool = True):
        self.use_lmfit = use_lmfit and LMFIT_AVAILABLE
        if not LMFIT_AVAILABLE and use_lmfit:
            logger.warning("lmfit not available, falling back to simple fitting")

    def fit_peak_advanced(self, x_data: np.ndarray, y_data: np.ndarray,
                          peak_center: float, model_type: str = 'voigt') -> Tuple[float, float, Dict]:
        """Advanced peak fitting using lmfit"""
        if not self.use_lmfit:
            return self.fit_peak_simple(x_data, y_data, peak_center)

        try:
            # Choose model based on type
            if model_type.lower() == 'voigt':
                model = models.VoigtModel()
            elif model_type.lower() == 'gaussian':
                model = models.GaussianModel()
            elif model_type.lower() == 'lorentzian':
                model = models.LorentzianModel()
            elif model_type.lower() == 'pseudo_voigt':
                model = models.PseudoVoigtModel()
            else:
                model = models.VoigtModel()  # Default

            # Add linear background
            model += models.LinearModel()

            # Set initial parameters
            params = model.make_params()
            params['center'].set(value=peak_center, min=peak_center - 0.5, max=peak_center + 0.5)
            params['amplitude'].set(value=np.max(y_data), min=0)

            if hasattr(params, 'sigma'):
                params['sigma'].set(value=0.1, min=0.01, max=1.0)
            if hasattr(params, 'gamma'):
                params['gamma'].set(value=0.1, min=0.01, max=1.0)

            # Fit
            result = model.fit(y_data, params, x=x_data)

            # Calculate FWHM based on model type
            fwhm = self._calculate_fwhm(result.params, model_type)

            # Calculate R-squared
            r_squared = 1 - result.residual.var() / np.var(y_data)

            # Additional fit information
            fit_info = {
                'model_type': model_type,
                'fit_params': {key: param.value for key, param in result.params.items()},
                'fit_errors': {key: param.stderr for key, param in result.params.items() if param.stderr is not None},
                'chi_squared': result.chisqr,
                'reduced_chi_squared': result.redchi,
                'aic': result.aic,
                'bic': result.bic,
                'fitted_curve': result.best_fit
            }

            return max(0.01, min(fwhm, 2.0)), max(0, min(r_squared, 1.0)), fit_info

        except Exception as e:
            logger.debug(f"Advanced peak fitting failed: {e}")
            return self.fit_peak_simple(x_data, y_data, peak_center)

    def _calculate_fwhm(self, params: Dict, model_type: str) -> float:
        """Calculate FWHM from fitted parameters"""
        try:
            if model_type.lower() == 'gaussian':
                sigma = params['sigma'].value
                return 2.355 * sigma  # 2*sqrt(2*ln(2)) * sigma
            elif model_type.lower() == 'lorentzian':
                gamma = params['gamma'].value
                return 2 * gamma
            elif model_type.lower() in ['voigt', 'pseudo_voigt']:
                if 'sigma' in params and 'gamma' in params:
                    sigma = params['sigma'].value
                    gamma = params['gamma'].value
                    # Voigt FWHM approximation
                    fg = 2 * sigma * np.sqrt(2 * np.log(2))  # Gaussian component
                    fl = 2 * gamma  # Lorentzian component
                    return 0.5346 * fl + np.sqrt(0.2166 * fl ** 2 + fg ** 2)
                elif 'sigma' in params:
                    return 2.355 * params['sigma'].value

            # Fallback
            return 0.1

        except Exception:
            return 0.1

    def fit_peak_simple(self, x_data: np.ndarray, y_data: np.ndarray,
                        peak_center: float) -> Tuple[float, float, Dict]:
        """Simple FWHM estimation using half-maximum method"""
        try:
            if len(x_data) < 5:
                return 0.1, 0.5, {'model_type': 'simple', 'method': 'half_maximum'}

            # Find peak maximum
            max_idx = np.argmax(y_data)
            peak_height = y_data[max_idx]

            # Estimate baseline
            edge_points = np.concatenate([y_data[:3], y_data[-3:]])
            baseline = np.median(edge_points) if len(edge_points) > 0 else np.min(y_data)

            half_max = baseline + (peak_height - baseline) / 2

            # Find half-maximum points
            left_idx = max_idx
            right_idx = max_idx

            while left_idx > 0 and y_data[left_idx] > half_max:
                left_idx -= 1
            while right_idx < len(y_data) - 1 and y_data[right_idx] > half_max:
                right_idx += 1

            if right_idx > left_idx:
                fwhm = x_data[right_idx] - x_data[left_idx]

                # Simple quality estimate based on peak shape
                peak_width_points = right_idx - left_idx
                if peak_width_points > 3:
                    r_squared = min(0.85, 0.5 + 0.1 * peak_width_points)
                else:
                    r_squared = 0.5

                fit_info = {
                    'model_type': 'simple',
                    'method': 'half_maximum',
                    'peak_height': peak_height,
                    'baseline': baseline,
                    'half_max': half_max,
                    'peak_width_points': peak_width_points
                }

                return max(0.05, min(fwhm, 2.0)), r_squared, fit_info
            else:
                return 0.1, 0.3, {'model_type': 'simple', 'method': 'failed'}

        except Exception as e:
            logger.debug(f"Simple peak fitting failed: {e}")
            return 0.1, 0.3, {'model_type': 'simple', 'method': 'error', 'error': str(e)}

    def fit_multiple_peaks(self, x_data: np.ndarray, y_data: np.ndarray,
                           peak_centers: list, model_type: str = 'voigt') -> Dict:
        """Fit multiple peaks simultaneously"""
        if not self.use_lmfit:
            # Fall back to individual peak fitting
            results = []
            for center in peak_centers:
                # Extract region around peak
                mask = np.abs(x_data - center) <= 1.0
                if np.sum(mask) >= 5:
                    fwhm, r2, info = self.fit_peak_simple(x_data[mask], y_data[mask], center)
                    results.append({'center': center, 'fwhm': fwhm, 'r_squared': r2, 'fit_info': info})
            return {'individual_fits': results, 'simultaneous_fit': False}

        try:
            # Build composite model
            composite_model = models.LinearModel()  # Background

            for i, center in enumerate(peak_centers):
                if model_type.lower() == 'voigt':
                    peak_model = models.VoigtModel(prefix=f'p{i}_')
                elif model_type.lower() == 'gaussian':
                    peak_model = models.GaussianModel(prefix=f'p{i}_')
                elif model_type.lower() == 'lorentzian':
                    peak_model = models.LorentzianModel(prefix=f'p{i}_')
                else:
                    peak_model = models.VoigtModel(prefix=f'p{i}_')

                composite_model += peak_model

            # Set initial parameters
            params = composite_model.make_params()

            for i, center in enumerate(peak_centers):
                params[f'p{i}_center'].set(value=center, min=center - 0.3, max=center + 0.3)
                params[f'p{i}_amplitude'].set(value=np.max(y_data) / len(peak_centers), min=0)

                if f'p{i}_sigma' in params:
                    params[f'p{i}_sigma'].set(value=0.1, min=0.01, max=0.5)
                if f'p{i}_gamma' in params:
                    params[f'p{i}_gamma'].set(value=0.1, min=0.01, max=0.5)

            # Fit
            result = composite_model.fit(y_data, params, x=x_data)

            # Extract individual peak results
            individual_results = []
            for i, center in enumerate(peak_centers):
                fwhm = self._calculate_fwhm_from_prefix(result.params, f'p{i}_', model_type)
                individual_results.append({
                    'center': result.params[f'p{i}_center'].value,
                    'fwhm': fwhm,
                    'amplitude': result.params[f'p{i}_amplitude'].value
                })

            r_squared = 1 - result.residual.var() / np.var(y_data)

            return {
                'individual_fits': individual_results,
                'simultaneous_fit': True,
                'overall_r_squared': r_squared,
                'chi_squared': result.chisqr,
                'fitted_curve': result.best_fit,
                'model_type': model_type
            }

        except Exception as e:
            logger.warning(f"Multiple peak fitting failed: {e}")
            # Fall back to individual fitting
            return self.fit_multiple_peaks(x_data, y_data, peak_centers, 'simple')

    def _calculate_fwhm_from_prefix(self, params: Dict, prefix: str, model_type: str) -> float:
        """Calculate FWHM from fitted parameters with prefix"""
        try:
            if model_type.lower() == 'gaussian':
                sigma = params[f'{prefix}sigma'].value
                return 2.355 * sigma
            elif model_type.lower() == 'lorentzian':
                gamma = params[f'{prefix}gamma'].value
                return 2 * gamma
            elif model_type.lower() == 'voigt':
                sigma = params[f'{prefix}sigma'].value
                gamma = params[f'{prefix}gamma'].value
                fg = 2 * sigma * np.sqrt(2 * np.log(2))
                fl = 2 * gamma
                return 0.5346 * fl + np.sqrt(0.2166 * fl ** 2 + fg ** 2)

            return 0.1

        except Exception:
            return 0.1