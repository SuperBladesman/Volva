"""
Standard Reference Material (SRM) module for XRD instrumental broadening correction.

This module provides comprehensive functionality for:
- Loading and processing SRM diffraction patterns
- Automatic SRM material identification
- Instrumental profile fitting with multiple functions
- Applying corrections to sample peak data
- Calibration management and validation

Usage:
    from volva.srm import SRMProcessor, InstrumentalProfileFitter, InstrumentalCorrector
    from volva.srm import get_srm, list_available_srms
"""

# Import main classes
from .processor import SRMProcessor, SRMProcessingResult, SRMPeakMatch
from .fitting import (
    InstrumentalProfileFitter,
    InstrumentalFitResult,
    CagliotiFunction,
    SimplePolynomialFunction,
    PearsonVIIFunction
)
from .correction import (
    InstrumentalCorrector,
    CorrectionCalibration,
    CorrectedPeak
)
from .standards import (
    StandardReferenceMaterial,
    get_srm,
    list_available_srms,
    recommend_srm_for_range
)

# Version info
__version__ = "1.0.0"

# Public API
__all__ = [
    # Main processing classes
    'SRMProcessor',
    'InstrumentalProfileFitter',
    'InstrumentalCorrector',

    # Data structures
    'SRMProcessingResult',
    'SRMPeakMatch',
    'InstrumentalFitResult',
    'CorrectionCalibration',
    'CorrectedPeak',
    'StandardReferenceMaterial',

    # Instrumental functions
    'CagliotiFunction',
    'SimplePolynomialFunction',
    'PearsonVIIFunction',

    # Utility functions
    'get_srm',
    'list_available_srms',
    'recommend_srm_for_range',

    # Convenience functions
    'create_srm_workflow',
    'quick_correction'
]


def create_srm_workflow(srm_file_path: str, wavelength: float,
                        srm_material: str = None) -> tuple:
    """
    Convenience function to create a complete SRM workflow.

    Args:
        srm_file_path: Path to SRM diffraction data
        wavelength: X-ray wavelength in Å
        srm_material: SRM material name (auto-detect if None)

    Returns:
        Tuple of (processor, corrector, calibration)
    """
    # Initialize components
    processor = SRMProcessor()
    corrector = InstrumentalCorrector()

    # Load and process SRM data
    two_theta, intensity, metadata = processor.load_srm_data(srm_file_path)

    # Identify material if not specified
    if srm_material is None:
        srm_material = processor.identify_srm_material(two_theta, intensity, wavelength)

    # Process SRM peaks
    srm_result = processor.process_srm_data(two_theta, intensity, srm_material, wavelength)

    # Create instrumental calibration
    instrument_config = metadata.get('instrument_config', {
        'source': 'Cu Ka',
        'voltage': 40,
        'current': 40
    })

    calibration = corrector.create_calibration(srm_result, instrument_config)

    return processor, corrector, calibration


def quick_correction(sample_peaks: list, calibration_file: str = None,
                     calibration: CorrectionCalibration = None) -> list:
    """
    Quick correction of sample peaks using existing calibration.

    Args:
        sample_peaks: List of sample peak fit results
        calibration_file: Path to saved calibration file
        calibration: Direct calibration object

    Returns:
        List of corrected peaks
    """
    corrector = InstrumentalCorrector()

    # Load calibration if file provided
    if calibration_file is not None:
        calibration = corrector.load_calibration(calibration_file)
    elif calibration is None:
        raise ValueError("Must provide either calibration_file or calibration object")

    # Apply corrections
    corrected_peaks = corrector.apply_correction(sample_peaks, calibration)

    return corrected_peaks


# Module-level convenience instances for simple usage
_default_processor = None
_default_corrector = None


def get_default_processor() -> SRMProcessor:
    """Get the default SRM processor instance."""
    global _default_processor
    if _default_processor is None:
        _default_processor = SRMProcessor()
    return _default_processor


def get_default_corrector() -> InstrumentalCorrector:
    """Get the default instrumental corrector instance."""
    global _default_corrector
    if _default_corrector is None:
        _default_corrector = InstrumentalCorrector()
    return _default_corrector