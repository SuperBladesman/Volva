"""
Volva XRD Analysis Tool - Exception Classes
Custom exception hierarchy for comprehensive error handling
"""


class VolvaError(Exception):
    """Base exception for Volva errors"""
    pass


class DataLoadError(VolvaError):
    """Exception raised when data loading fails"""
    pass


class PeakDetectionError(VolvaError):
    """Exception raised when peak detection fails"""
    pass


class AnalysisError(VolvaError):
    """Exception raised when analysis fails"""
    pass


class InsufficientDataError(AnalysisError):
    """Exception raised when insufficient data is provided for analysis"""
    pass


class ConfigurationError(VolvaError):
    """Exception raised when configuration is invalid"""
    pass


class ValidationError(VolvaError):
    """Exception raised when data validation fails"""
    pass


class FittingError(AnalysisError):
    """Exception raised when peak fitting fails"""
    pass


class ExportError(VolvaError):
    """Exception raised when data export fails"""
    pass

# Aliases for backward compatibility and consistency
VolvaDataError = DataLoadError
VolvaPeakError = PeakDetectionError
VolvaFittingError = FittingError
VolvaAnalysisError = AnalysisError
VolvaConfigError = ConfigurationError
VolvaValidationError = ValidationError
VolvaExportError = ExportError
VolvaCorrectionError = AnalysisError