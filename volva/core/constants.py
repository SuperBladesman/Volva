"""
Volva XRD Analysis Tool - Physical Constants and Wavelength Data
Contains X-ray wavelengths, Scherrer constants, and other physical constants
"""

import logging
from typing import Dict

logger = logging.getLogger(__name__)


class Constants:
    """Physical constants and wavelength data with extensibility"""

    # Standard wavelengths - can be extended
    WAVELENGTHS = {
        "Cu Ka": {"ka1": 0.15406, "ka2": 0.154439, "ratio": 0.5, "alpha": 2.2e-3},
        "Mo Ka": {"ka1": 0.07093, "ka2": 0.071073, "ratio": 0.5, "alpha": 1.86e-3},
        "Ag Ka": {"ka1": 0.05594, "ka2": 0.05638, "ratio": 0.5, "alpha": 1.54e-3},
        "Co Ka": {"ka1": 0.178897, "ka2": 0.179285, "ratio": 0.5, "alpha": 2.39e-3},
        "Cr Ka": {"ka1": 0.228970, "ka2": 0.229361, "ratio": 0.5, "alpha": 2.93e-3},
        "Fe Ka": {"ka1": 0.193604, "ka2": 0.193998, "ratio": 0.5, "alpha": 2.56e-3},
        "W La": {"ka1": 0.147640, "ka2": 0.148480, "ratio": 0.5, "alpha": 2.1e-3},
    }

    # Scherrer constants with better documentation
    SCHERRER_CONSTANTS = {
        "General/Approximation": {"value": 0.9, "description": "General approximation"},
        "Spherical (FWHM)": {"value": 0.94, "description": "Spherical particles, FWHM"},
        "Spherical (Integral Breadth)": {"value": 0.89, "description": "Spherical particles, integral breadth"},
        "Cubic": {"value": 1.0, "description": "Cubic crystallites"},
        "Rod-like": {"value": 1.84, "description": "Rod-like crystallites"},
        "Plate-like": {"value": 0.89, "description": "Plate-like crystallites"},
        "Tetragonal": {"value": 0.89, "description": "Tetragonal crystallites"},
    }

    @classmethod
    def add_wavelength(cls, name: str, ka1: float, ka2: float = None,
                       ratio: float = 0.5, alpha: float = 2.0e-3):
        """Add new wavelength to the database"""
        if ka2 is None:
            ka2 = ka1 * 1.002  # Approximate Ka2 from Ka1

        cls.WAVELENGTHS[name] = {
            "ka1": ka1, "ka2": ka2, "ratio": ratio, "alpha": alpha
        }
        logger.info(f"Added wavelength: {name} (λ = {ka1:.5f} nm)")

    @classmethod
    def add_scherrer_constant(cls, name: str, value: float, description: str = ""):
        """Add new Scherrer constant"""
        cls.SCHERRER_CONSTANTS[name] = {"value": value, "description": description}
        logger.info(f"Added Scherrer constant: {name} = {value}")

    @classmethod
    def get_wavelength_info(cls, name: str) -> Dict:
        """Get wavelength information by name"""
        if name not in cls.WAVELENGTHS:
            raise ValueError(f"Unknown wavelength: {name}")
        return cls.WAVELENGTHS[name].copy()

    @classmethod
    def get_scherrer_constant(cls, name: str) -> Dict:
        """Get Scherrer constant information by name"""
        if name not in cls.SCHERRER_CONSTANTS:
            raise ValueError(f"Unknown Scherrer constant: {name}")
        return cls.SCHERRER_CONSTANTS[name].copy()

    @classmethod
    def list_wavelengths(cls) -> list:
        """List all available wavelength names"""
        return list(cls.WAVELENGTHS.keys())

    @classmethod
    def list_scherrer_constants(cls) -> list:
        """List all available Scherrer constant names"""
        return list(cls.SCHERRER_CONSTANTS.keys())

# Export class attributes as module-level variables for backward compatibility
WAVELENGTHS = Constants.WAVELENGTHS
SCHERRER_CONSTANTS = Constants.SCHERRER_CONSTANTS