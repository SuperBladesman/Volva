"""
Volva XRD Analysis Tool - Core Package
Contains fundamental classes: configuration, constants, and exceptions
"""

import logging

# Configure logging for the core package
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Import main classes for easy access
from .constants import Constants
from .config import VolvaConfig
from .exceptions import (
    VolvaError,
    DataLoadError,
    PeakDetectionError,
    AnalysisError,
    InsufficientDataError,
    ConfigurationError,
    ValidationError,
    FittingError,
    ExportError
)

# Check for optional dependencies
try:
    import lmfit
    LMFIT_AVAILABLE = True
except ImportError:
    LMFIT_AVAILABLE = False
    logging.getLogger(__name__).warning("lmfit not available. Using simplified peak fitting.")

# Package version
__version__ = "4.0.0"

# Public API
__all__ = [
    'Constants',
    'VolvaConfig',
    'VolvaError',
    'DataLoadError',
    'PeakDetectionError',
    'AnalysisError',
    'InsufficientDataError',
    'ConfigurationError',
    'ValidationError',
    'FittingError',
    'ExportError',
    'LMFIT_AVAILABLE',
    '__version__'
]