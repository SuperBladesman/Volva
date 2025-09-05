"""
Instrumental profile fitting for SRM-based broadening correction.

Implements various instrumental broadening functions and fitting procedures
to determine angle-dependent instrumental contribution to peak widths.
"""

import numpy as np
from scipy.optimize import curve_fit, minimize
from typing import Dict, List, Tuple, Optional, Callable, NamedTuple
import warnings
from abc import ABC, abstractmethod

from .processor import SRMProcessingResult, SRMPeakMatch
from ..core.exceptions import VolvaFittingError


class InstrumentalFitResult(NamedTuple):
    """Results from instrumental profile fitting."""
    function_name: str
    parameters: Dict[str, float]
    parameter_errors: Dict[str, float]
    r_squared: float
    chi_squared: float
    residuals: np.ndarray
    fitted_values: np.ndarray
    two_theta_range: Tuple[float, float]


class InstrumentalFunction(ABC):
    """Base class for instrumental broadening functions."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Function name."""
        pass

    @property
    @abstractmethod
    def parameter_names(self) -> List[str]:
        """List of parameter names."""
        pass

    @property
    @abstractmethod
    def initial_guess(self) -> List[float]:
        """Initial parameter guess."""
        pass

    @abstractmethod
    def function(self, two_theta: np.ndarray, *params) -> np.ndarray:
        """
        Calculate instrumental FWHM as function of 2θ.

        Args:
            two_theta: 2θ angles in degrees
            *params: Function parameters

        Returns:
            FWHM values in degrees
        """
        pass

    @abstractmethod
    def bounds(self) -> Tuple[List[float], List[float]]:
        """Parameter bounds as (lower, upper) lists."""
        pass


class CagliotiFunction(InstrumentalFunction):
    """
    Caglioti function: FWHM² = U·tan²θ + V·tanθ + W

    Most commonly used for conventional powder diffractometers.
    U, V, W are the Caglioti parameters.
    """

    @property
    def name(self) -> str:
        return "Caglioti"

    @property
    def parameter_names(self) -> List[str]:
        return ["U", "V", "W"]

    @property
    def initial_guess(self) -> List[float]:
        return [0.01, -0.005, 0.005]  # Typical starting values

    def function(self, two_theta: np.ndarray, U: float, V: float, W: float) -> np.ndarray:
        """Calculate FWHM using Caglioti function."""
        theta_rad = np.radians(two_theta / 2)
        tan_theta = np.tan(theta_rad)

        fwhm_squared = U * tan_theta ** 2 + V * tan_theta + W

        # Ensure positive FWHM values
        fwhm_squared = np.maximum(fwhm_squared, 1e-6)

        return np.sqrt(fwhm_squared)

    def bounds(self) -> Tuple[List[float], List[float]]:
        return ([-0.1, -0.1, 0.0], [0.1, 0.1, 0.1])


class SimplePolynomialFunction(InstrumentalFunction):
    """
    Simple polynomial: FWHM = A + B·2θ + C·(2θ)²

    Alternative to Caglioti for some instruments.
    """

    @property
    def name(self) -> str:
        return "Polynomial"

    @property
    def parameter_names(self) -> List[str]:
        return ["A", "B", "C"]

    @property
    def initial_guess(self) -> List[float]:
        return [0.1, 0.001, 0.00001]

    def function(self, two_theta: np.ndarray, A: float, B: float, C: float) -> np.ndarray:
        """Calculate FWHM using polynomial function."""
        fwhm = A + B * two_theta + C * two_theta ** 2
        return np.maximum(fwhm, 0.001)  # Minimum FWHM to avoid negatives

    def bounds(self) -> Tuple[List[float], List[float]]:
        return ([0.001, -0.01, -0.001], [1.0, 0.01, 0.001])


class PearsonVIIFunction(InstrumentalFunction):
    """
    Modified Pearson VII approach: FWHM = A/(cos²θ) + B·tanθ + C

    Good for synchrotron and some laboratory sources.
    """

    @property
    def name(self) -> str:
        return "Pearson-VII"

    @property
    def parameter_names(self) -> List[str]:
        return ["A", "B", "C"]

    @property
    def initial_guess(self) -> List[float]:
        return [0.001, 0.005, 0.01]

    def function(self, two_theta: np.ndarray, A: float, B: float, C: float) -> np.ndarray:
        """Calculate FWHM using modified Pearson VII function."""
        theta_rad = np.radians(two_theta / 2)
        cos_theta = np.cos(theta_rad)
        tan_theta = np.tan(theta_rad)

        # Avoid division by zero
        cos_theta = np.maximum(cos_theta, 1e-6)

        fwhm = A / cos_theta ** 2 + B * tan_theta + C
        return np.maximum(fwhm, 0.001)

    def bounds(self) -> Tuple[List[float], List[float]]:
        return ([0.0, 0.0, 0.001], [0.1, 0.1, 0.5])


class InstrumentalProfileFitter:
    """Fits instrumental broadening functions to SRM peak width data."""

    def __init__(self):
        """Initialize the fitter with available functions."""
        self.functions = {
            'Caglioti': CagliotiFunction(),
            'Polynomial': SimplePolynomialFunction(),
            'Pearson-VII': PearsonVIIFunction()
        }
        self.default_function = 'Caglioti'

    def fit_instrumental_profile(self, srm_result: SRMProcessingResult,
                                 function_name: Optional[str] = None,
                                 weight_by_intensity: bool = True) -> InstrumentalFitResult:
        """
        Fit instrumental broadening function to SRM peak widths.

        Args:
            srm_result: Processed SRM data with matched peaks
            function_name: Name of function to fit (default: Caglioti)
            weight_by_intensity: Whether to weight fit by peak intensities

        Returns:
            Instrumental fit results
        """
        if len(srm_result.matched_peaks) < 3:
            raise VolvaFittingError("Insufficient peaks for instrumental fitting (need ≥3)")

        function_name = function_name or self.default_function
        if function_name not in self.functions:
            raise VolvaFittingError(f"Unknown function: {function_name}")

        func_obj = self.functions[function_name]

        # Extract peak data
        two_theta_data = np.array([peak.observed_position for peak in srm_result.matched_peaks])
        fwhm_data = np.array([peak.fwhm for peak in srm_result.matched_peaks])

        # Validate FWHM data
        if np.any(fwhm_data <= 0):
            raise VolvaFittingError("Invalid FWHM values (≤ 0) in SRM data")

        # Prepare weights
        if weight_by_intensity:
            intensities = np.array([peak.intensity for peak in srm_result.matched_peaks])
            weights = np.sqrt(intensities) / np.max(intensities)  # Normalize weights
        else:
            weights = np.ones_like(fwhm_data)

        # Perform fitting
        try:
            popt, pcov = curve_fit(
                func_obj.function,
                two_theta_data,
                fwhm_data,
                p0=func_obj.initial_guess,
                bounds=func_obj.bounds(),
                sigma=1.0 / weights,
                absolute_sigma=False,
                maxfev=5000
            )
        except Exception as e:
            raise VolvaFittingError(f"Fitting failed: {str(e)}")

        # Calculate fit statistics
        fitted_fwhm = func_obj.function(two_theta_data, *popt)
        residuals = fwhm_data - fitted_fwhm

        # R-squared
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((fwhm_data - np.mean(fwhm_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        # Chi-squared
        chi_squared = np.sum((residuals / weights) ** 2) / (len(fwhm_data) - len(popt))

        # Parameter errors
        try:
            param_errors = np.sqrt(np.diag(pcov))
        except:
            param_errors = np.full_like(popt, np.nan)

        # Build results
        parameters = dict(zip(func_obj.parameter_names, popt))
        parameter_errors = dict(zip(func_obj.parameter_names, param_errors))

        return InstrumentalFitResult(
            function_name=function_name,
            parameters=parameters,
            parameter_errors=parameter_errors,
            r_squared=r_squared,
            chi_squared=chi_squared,
            residuals=residuals,
            fitted_values=fitted_fwhm,
            two_theta_range=(np.min(two_theta_data), np.max(two_theta_data))
        )

    def compare_functions(self, srm_result: SRMProcessingResult,
                          functions: Optional[List[str]] = None) -> Dict[str, InstrumentalFitResult]:
        """
        Compare multiple instrumental functions and rank by fit quality.

        Args:
            srm_result: Processed SRM data
            functions: List of function names to compare (all if None)

        Returns:
            Dictionary of fit results for each function
        """
        if functions is None:
            functions = list(self.functions.keys())

        results = {}

        for func_name in functions:
            try:
                result = self.fit_instrumental_profile(srm_result, func_name)
                results[func_name] = result
            except VolvaFittingError as e:
                print(f"Warning: Failed to fit {func_name}: {str(e)}")
                continue

        if not results:
            raise VolvaFittingError("No functions could be fitted successfully")

        return results

    def select_best_function(self, fit_results: Dict[str, InstrumentalFitResult],
                             criterion: str = 'r_squared') -> str:
        """
        Select best function based on fit quality criterion.

        Args:
            fit_results: Dictionary of fit results from compare_functions
            criterion: Selection criterion ('r_squared', 'chi_squared', or 'aic')

        Returns:
            Name of best function
        """
        if not fit_results:
            raise VolvaFittingError("No fit results provided")

        if criterion == 'r_squared':
            # Higher R² is better
            best_func = max(fit_results.keys(),
                            key=lambda k: fit_results[k].r_squared)
        elif criterion == 'chi_squared':
            # Lower χ² is better
            best_func = min(fit_results.keys(),
                            key=lambda k: fit_results[k].chi_squared)
        elif criterion == 'aic':
            # Akaike Information Criterion (lower is better)
            def aic_score(result):
                n = len(result.fitted_values)
                k = len(result.parameters)
                return 2 * k + n * np.log(np.sum(result.residuals ** 2))

            best_func = min(fit_results.keys(),
                            key=lambda k: aic_score(fit_results[k]))
        else:
            raise ValueError(f"Unknown criterion: {criterion}")

        return best_func

    def predict_instrumental_fwhm(self, two_theta: np.ndarray,
                                  fit_result: InstrumentalFitResult) -> np.ndarray:
        """
        Predict instrumental FWHM at given 2θ positions.

        Args:
            two_theta: 2θ positions in degrees
            fit_result: Results from instrumental fitting

        Returns:
            Predicted FWHM values in degrees
        """
        func_obj = self.functions[fit_result.function_name]
        params = [fit_result.parameters[name] for name in func_obj.parameter_names]

        return func_obj.function(two_theta, *params)

    def validate_fit_quality(self, fit_result: InstrumentalFitResult,
                             min_r_squared: float = 0.8,
                             max_chi_squared: float = 5.0) -> Dict[str, bool]:
        """
        Validate the quality of an instrumental fit.

        Args:
            fit_result: Fitting results to validate
            min_r_squared: Minimum acceptable R²
            max_chi_squared: Maximum acceptable χ²

        Returns:
            Dictionary of validation results
        """
        validation = {
            'r_squared_ok': fit_result.r_squared >= min_r_squared,
            'chi_squared_ok': fit_result.chi_squared <= max_chi_squared,
            'no_negative_fwhm': np.all(fit_result.fitted_values > 0),
            'reasonable_residuals': np.std(fit_result.residuals) < 0.05,  # < 0.05° std
            'parameters_reasonable': self._check_parameter_reasonableness(fit_result)
        }

        validation['overall_ok'] = all(validation.values())

        return validation

    def _check_parameter_reasonableness(self, fit_result: InstrumentalFitResult) -> bool:
        """Check if fitted parameters are within reasonable ranges."""
        func_name = fit_result.function_name
        params = fit_result.parameters

        if func_name == 'Caglioti':
            # Typical ranges for Caglioti parameters
            return (abs(params['U']) < 0.1 and
                    abs(params['V']) < 0.1 and
                    params['W'] > 0 and params['W'] < 0.1)
        elif func_name == 'Polynomial':
            # Check for reasonable polynomial coefficients
            return (params['A'] > 0 and params['A'] < 1.0 and
                    abs(params['B']) < 0.01 and
                    abs(params['C']) < 0.001)
        elif func_name == 'Pearson-VII':
            # Check Pearson VII parameters
            return (params['A'] >= 0 and params['A'] < 0.1 and
                    params['B'] >= 0 and params['B'] < 0.1 and
                    params['C'] > 0 and params['C'] < 0.5)

        return True  # Default to reasonable if unknown function