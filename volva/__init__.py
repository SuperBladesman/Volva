"""
Volva XRD Analysis Tool v4.0 Enhanced
Advanced X-ray Diffraction analysis with modular architecture and scientific validation
"""

# Import core functionality
from .core import (
    VolvaConfig, Constants, LMFIT_AVAILABLE,
    VolvaError, DataLoadError, PeakDetectionError, AnalysisError,
    InsufficientDataError, ConfigurationError, ValidationError,
    FittingError, ExportError
)

# Import data processing
from .data import XRDDataProcessor

# Import peak detection and fitting
from .peaks import PeakModel, PeakFitter, PeakDetector

# Import analysis methods
from .analysis import (
    AnalysisMethod, StrainAnalysisMethod,
    ScherrerAnalyzer, MonshiScherrerAnalyzer,
    WilliamsonHallAnalyzer, HalderWagnerAnalyzer
)

# Import interfaces
from .interfaces import JupyterVolva, VolvaGUI, launch_gui

# Version information
__version__ = "4.0.0"
__author__ = "Volva Development Team"
__license__ = "MIT"

# Public API
__all__ = [
    # Version info
    '__version__',

    # Core components
    'VolvaConfig',
    'Constants',
    'LMFIT_AVAILABLE',

    # Exceptions
    'VolvaError',
    'DataLoadError',
    'PeakDetectionError',
    'AnalysisError',
    'InsufficientDataError',
    'ConfigurationError',
    'ValidationError',
    'FittingError',
    'ExportError',

    # Data processing
    'XRDDataProcessor',

    # Peak analysis
    'PeakModel',
    'PeakFitter',
    'PeakDetector',

    # Analysis methods
    'AnalysisMethod',
    'StrainAnalysisMethod',
    'ScherrerAnalyzer',
    'MonshiScherrerAnalyzer',
    'WilliamsonHallAnalyzer',
    'HalderWagnerAnalyzer',

    # Interfaces
    'JupyterVolva',
    'VolvaGUI',
    'launch_gui'
]


# Convenience functions
def create_analyzer(wavelength: str = "Cu Ka", config_file: str = None):
    """
    Create a Jupyter-friendly analyzer instance

    Parameters:
    - wavelength: X-ray source name (e.g., "Cu Ka", "Mo Ka")
    - config_file: Optional configuration file path

    Returns:
    - JupyterVolva instance
    """
    return JupyterVolva(wavelength, config_file)


def quick_analysis(file_path: str, wavelength: str = "Cu Ka", config_file: str = None):
    """
    Perform quick comprehensive analysis

    Parameters:
    - file_path: Path to XRD data file
    - wavelength: X-ray source name
    - config_file: Optional configuration file path

    Returns:
    - Tuple of (analyzer, results) or (analyzer, None) if failed
    """
    try:
        analyzer = create_analyzer(wavelength, config_file)

        if analyzer.load_data(file_path):
            peaks = analyzer.detect_peaks()
            if peaks:
                results = analyzer.run_all_methods(verbose=False)
                return analyzer, results

        return analyzer, None

    except Exception as e:
        print(f"Quick analysis failed: {e}")
        return None, None


def show_info():
    """Display Volva information"""
    info_text = f"""
🔬 Volva XRD Analysis Tool v{__version__}
{'=' * 50}

📊 ANALYSIS METHODS:
• Scherrer Analysis - Basic crystallite size
• Monshi-Scherrer - Enhanced with Kα2 correction
• Williamson-Hall - Size + microstrain (isotropic)
• Halder-Wagner - Size + microstrain (anisotropic)

✨ ENHANCED FEATURES:
• Modular architecture with robust error handling
• Advanced peak detection with lmfit integration
• Scientific validation of analysis assumptions
• Quality-based filtering and outlier detection
• Comprehensive uncertainty quantification
• Multiple export formats

🚀 QUICK START:
# Jupyter Interface
import volva
analyzer = volva.create_analyzer('Cu Ka')
analyzer.load_data('your_file.txt')
results = analyzer.run_all_methods()

# GUI Interface  
volva.launch_gui()

# Quick Analysis
analyzer, results = volva.quick_analysis('your_file.txt')

🧮 Advanced Fitting: {'Available (lmfit)' if LMFIT_AVAILABLE else 'Basic only (install lmfit)'}
📡 Available Sources: {len(Constants.WAVELENGTHS)} X-ray wavelengths
⚙️  Shape Factors: {len(Constants.SCHERRER_CONSTANTS)} Scherrer constants

For detailed help: help(volva.JupyterVolva)
For GUI interface: volva.launch_gui()
    """
    print(info_text)


# Display info on import
print(f"🔬 Volva XRD Analysis Tool v{__version__} - Enhanced Edition")
print(f"📊 Ready for crystallite size and microstrain analysis")
print(f"🧮 Advanced fitting: {'Available' if LMFIT_AVAILABLE else 'Basic (install lmfit for more features)'}")
print(f"💡 Use volva.show_info() for detailed information")