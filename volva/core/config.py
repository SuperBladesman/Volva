"""
Volva XRD Analysis Tool - Configuration Management
Enhanced configuration class for Volva application parameters
"""

import json
import logging
from dataclasses import dataclass, field
from typing import Dict

from .exceptions import ConfigurationError

logger = logging.getLogger(__name__)


@dataclass
class VolvaConfig:
    """Enhanced configuration class for Volva application"""

    # Peak detection parameters
    proximity_tolerance_deg: float = 5.0
    default_bg_filter_window: int = 151
    default_prominence_factor: float = 0.1
    default_min_width_deg: float = 0.05
    default_fitting_window_deg: float = 1.0

    # Analysis parameters
    default_scherrer_k: str = "Spherical (FWHM)"
    default_wh_k: float = 0.94  # Made explicit and configurable
    default_hw_k: float = 0.94  # Made explicit and configurable
    default_hw_a: float = 2.0
    default_wavelength: str = "Cu Ka"
    min_r_squared: float = 0.8

    # Background subtraction
    default_bg_method: str = "polynomial"
    default_poly_order: int = 3

    # Peak fitting
    default_peak_model: str = "voigt"
    use_lmfit: bool = True

    # Instrumental broadening (configurable per instrument)
    default_instrumental_broadening: float = 0.05
    instrumental_broadening_models: Dict[str, float] = field(default_factory=lambda: {
        "Laboratory Diffractometer": 0.05,
        "Synchrotron": 0.01,
        "Neutron": 0.03,
        "Portable XRD": 0.1,
        "Custom": 0.05
    })

    def to_dict(self) -> Dict:
        """Convert config to dictionary"""
        return {
            'proximity_tolerance_deg': self.proximity_tolerance_deg,
            'default_bg_filter_window': self.default_bg_filter_window,
            'default_prominence_factor': self.default_prominence_factor,
            'default_min_width_deg': self.default_min_width_deg,
            'default_fitting_window_deg': self.default_fitting_window_deg,
            'default_scherrer_k': self.default_scherrer_k,
            'default_wh_k': self.default_wh_k,
            'default_hw_k': self.default_hw_k,
            'default_hw_a': self.default_hw_a,
            'default_wavelength': self.default_wavelength,
            'min_r_squared': self.min_r_squared,
            'default_bg_method': self.default_bg_method,
            'default_poly_order': self.default_poly_order,
            'default_peak_model': self.default_peak_model,
            'use_lmfit': self.use_lmfit,
            'default_instrumental_broadening': self.default_instrumental_broadening,
            'instrumental_broadening_models': self.instrumental_broadening_models
        }

    @classmethod
    def from_dict(cls, data: Dict) -> 'VolvaConfig':
        """Create config from dictionary"""
        try:
            return cls(**data)
        except TypeError as e:
            raise ConfigurationError(f"Invalid configuration data: {e}")

    def save(self, filepath: str):
        """Save configuration to file"""
        try:
            with open(filepath, 'w') as f:
                json.dump(self.to_dict(), f, indent=2)
            logger.info(f"Configuration saved to {filepath}")
        except Exception as e:
            logger.error(f"Failed to save configuration: {e}")
            raise ConfigurationError(f"Configuration save failed: {e}")

    @classmethod
    def load(cls, filepath: str) -> 'VolvaConfig':
        """Load configuration from file"""
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            logger.info(f"Configuration loaded from {filepath}")
            return cls.from_dict(data)
        except FileNotFoundError:
            logger.warning(f"Configuration file {filepath} not found, using defaults")
            return cls()
        except Exception as e:
            logger.error(f"Failed to load configuration: {e}")
            raise ConfigurationError(f"Configuration load failed: {e}")

    def validate(self) -> bool:
        """Validate configuration parameters"""
        errors = []

        # Validate numeric ranges
        if not (0.001 <= self.default_prominence_factor <= 1.0):
            errors.append("default_prominence_factor must be between 0.001 and 1.0")

        if not (0.01 <= self.default_min_width_deg <= 5.0):
            errors.append("default_min_width_deg must be between 0.01 and 5.0")

        if not (0.5 <= self.min_r_squared <= 1.0):
            errors.append("min_r_squared must be between 0.5 and 1.0")

        if not (0.001 <= self.default_instrumental_broadening <= 1.0):
            errors.append("default_instrumental_broadening must be between 0.001 and 1.0")

        # Validate string parameters
        valid_bg_methods = ["polynomial", "rolling_minimum", "linear"]
        if self.default_bg_method not in valid_bg_methods:
            errors.append(f"default_bg_method must be one of {valid_bg_methods}")

        valid_peak_models = ["gaussian", "lorentzian", "voigt"]
        if self.default_peak_model not in valid_peak_models:
            errors.append(f"default_peak_model must be one of {valid_peak_models}")

        if errors:
            raise ConfigurationError("Configuration validation failed: " + "; ".join(errors))

        return True

    def update_instrumental_broadening(self, instrument_type: str):
        """Update instrumental broadening based on instrument type"""
        if instrument_type in self.instrumental_broadening_models:
            self.default_instrumental_broadening = self.instrumental_broadening_models[instrument_type]
            logger.info(
                f"Updated instrumental broadening to {self.default_instrumental_broadening:.4f}° for {instrument_type}")
        else:
            logger.warning(f"Unknown instrument type: {instrument_type}")

    def add_instrument_model(self, name: str, broadening: float):
        """Add new instrument model"""
        self.instrumental_broadening_models[name] = broadening
        logger.info(f"Added instrument model: {name} with broadening {broadening:.4f}°")