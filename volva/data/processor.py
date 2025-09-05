"""
Volva XRD Analysis Tool - Data Processing Module
Enhanced XRD data loading and processing with comprehensive error handling
"""

import os
import numpy as np
from scipy.signal import savgol_filter
from scipy.ndimage import minimum_filter1d
from typing import Tuple, Optional
import logging

from ..core import DataLoadError

logger = logging.getLogger(__name__)


class XRDDataProcessor:
    """Enhanced XRD data loading and processing with better error handling"""

    @staticmethod
    def load_xrd_file(file_path: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], str]:
        """Load XRD data with comprehensive format detection and validation"""
        if not os.path.exists(file_path):
            raise DataLoadError(f"File not found: {file_path}")

        if not os.path.isfile(file_path):
            raise DataLoadError(f"Path is not a file: {file_path}")

        try:
            # Try different delimiters and skip configurations
            delimiters = [None, '\t', ',', ' ', ';']
            skip_rows_options = [0, 1, 2, 3]

            for delimiter in delimiters:
                for skip_rows in skip_rows_options:
                    try:
                        data = np.loadtxt(file_path, delimiter=delimiter, skiprows=skip_rows)

                        if data.ndim == 1:
                            # Single column - assume intensity only
                            if len(data) < 10:
                                continue
                            two_theta = np.linspace(10, 80, len(data))
                            intensity = data
                        elif data.ndim == 2:
                            if data.shape[1] < 2:
                                continue
                            two_theta = data[:, 0]
                            intensity = data[:, 1]
                        else:
                            continue

                        # Validate data
                        if len(two_theta) < 10:
                            continue

                        if not (5 <= np.min(two_theta) <= 15 and 70 <= np.max(two_theta) <= 180):
                            # Try swapping columns if range is wrong
                            if data.ndim == 2 and data.shape[1] >= 2:
                                two_theta, intensity = intensity, two_theta
                                if not (5 <= np.min(two_theta) <= 15 and 70 <= np.max(two_theta) <= 180):
                                    continue

                        # Final validation
                        if np.any(np.isnan(two_theta)) or np.any(np.isnan(intensity)):
                            continue

                        if np.any(intensity < 0):
                            intensity = np.maximum(intensity, 0)  # Remove negative values

                        # Sort by 2theta
                        sort_indices = np.argsort(two_theta)
                        two_theta = two_theta[sort_indices]
                        intensity = intensity[sort_indices]

                        # Remove duplicates
                        unique_indices = np.unique(two_theta, return_index=True)[1]
                        two_theta = two_theta[unique_indices]
                        intensity = intensity[unique_indices]

                        success_msg = f"Data loaded successfully: {len(two_theta)} points"
                        success_msg += f" (Range: {two_theta.min():.1f}° - {two_theta.max():.1f}°)"
                        success_msg += f" (Skip rows: {skip_rows}, Delimiter: {delimiter or 'whitespace'})"

                        logger.info(success_msg)
                        return two_theta, intensity, success_msg

                    except (ValueError, IndexError, UnicodeDecodeError):
                        continue

            raise DataLoadError("Could not parse data file with any known format")

        except Exception as e:
            if isinstance(e, DataLoadError):
                raise
            raise DataLoadError(f"Error loading file: {str(e)}")

    @staticmethod
    def subtract_background(two_theta: np.ndarray, intensity: np.ndarray,
                            method: str = 'polynomial', poly_order: int = 3,
                            **kwargs) -> Tuple[np.ndarray, np.ndarray]:
        """Enhanced background subtraction with multiple methods"""
        try:
            if method == 'polynomial':
                # Polynomial background with outlier rejection
                coeffs = np.polyfit(two_theta, intensity, poly_order)
                background = np.polyval(coeffs, two_theta)

                # Iterative outlier rejection
                for _ in range(3):
                    residuals = intensity - background
                    threshold = np.percentile(np.abs(residuals), 85)
                    mask = np.abs(residuals) < threshold
                    if np.sum(mask) < len(intensity) * 0.5:
                        break
                    coeffs = np.polyfit(two_theta[mask], intensity[mask], poly_order)
                    background = np.polyval(coeffs, two_theta)

            elif method == 'rolling_minimum':
                # Rolling minimum with smoothing
                window = kwargs.get('window', min(151, len(intensity) // 5))
                if window % 2 == 0:
                    window += 1
                background = minimum_filter1d(intensity, size=window, mode='nearest')

                # Smooth the background
                if len(background) > 20:
                    smooth_window = min(51, max(5, len(background) // 20))
                    if smooth_window % 2 == 0:
                        smooth_window += 1
                    background = savgol_filter(background, smooth_window, 2)

            elif method == 'linear':
                # Linear background between first and last points
                slope = (intensity[-1] - intensity[0]) / (two_theta[-1] - two_theta[0])
                background = intensity[0] + slope * (two_theta - two_theta[0])

            else:
                raise ValueError(f"Unknown background method: {method}")

            # Ensure background doesn't exceed intensity
            background = np.minimum(background, intensity * 0.98)
            background_subtracted = intensity - background
            background_subtracted = np.maximum(background_subtracted, 0)

            logger.info(f"Background subtraction completed using {method} method")
            return background_subtracted, background

        except Exception as e:
            logger.warning(f"Background subtraction failed: {e}")
            return intensity.copy(), np.zeros_like(intensity)

    @staticmethod
    def validate_xrd_data(two_theta: np.ndarray, intensity: np.ndarray) -> bool:
        """Validate XRD data quality and range"""
        try:
            # Check for valid arrays
            if two_theta is None or intensity is None:
                return False

            if len(two_theta) != len(intensity):
                return False

            if len(two_theta) < 10:
                return False

            # Check for reasonable 2theta range
            if not (5 <= np.min(two_theta) <= 20 and 60 <= np.max(two_theta) <= 180):
                return False

            # Check for NaN or infinite values
            if np.any(~np.isfinite(two_theta)) or np.any(~np.isfinite(intensity)):
                return False

            # Check for negative intensities (shouldn't be many)
            if np.sum(intensity < 0) > len(intensity) * 0.1:
                return False

            return True

        except Exception as e:
            logger.warning(f"Data validation failed: {e}")
            return False

    @staticmethod
    def smooth_data(intensity: np.ndarray, window_length: int = None, poly_order: int = 3) -> np.ndarray:
        """Smart data smoothing with Savitzky-Golay filter"""
        try:
            if len(intensity) <= 10:
                return intensity.copy()

            if window_length is None:
                window_length = min(51, max(5, len(intensity) // 20))

            # Ensure odd window length
            if window_length % 2 == 0:
                window_length += 1

            # Ensure window length is not larger than data
            window_length = min(window_length, len(intensity))
            if window_length < poly_order + 1:
                window_length = poly_order + 1
                if window_length % 2 == 0:
                    window_length += 1

            return savgol_filter(intensity, window_length, poly_order)

        except Exception as e:
            logger.warning(f"Data smoothing failed: {e}")
            return intensity.copy()

    @staticmethod
    def resample_data(two_theta: np.ndarray, intensity: np.ndarray,
                      target_points: int = None, target_step: float = None) -> Tuple[np.ndarray, np.ndarray]:
        """Resample XRD data to uniform spacing"""
        try:
            if target_points is None and target_step is None:
                target_points = len(two_theta)

            if target_step is not None:
                # Calculate points needed for given step size
                range_2theta = two_theta.max() - two_theta.min()
                target_points = int(range_2theta / target_step) + 1

            # Create uniform 2theta array
            two_theta_uniform = np.linspace(two_theta.min(), two_theta.max(), target_points)

            # Interpolate intensity values
            intensity_uniform = np.interp(two_theta_uniform, two_theta, intensity)

            logger.info(f"Data resampled from {len(two_theta)} to {len(two_theta_uniform)} points")
            return two_theta_uniform, intensity_uniform

        except Exception as e:
            logger.warning(f"Data resampling failed: {e}")
            return two_theta, intensity