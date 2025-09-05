"""
Volva XRD Analysis Tool - Peaks Package
Peak detection, fitting, and analysis functionality
"""

from .models import PeakModel, PeakFitter, PeakFitResult
from .detector import PeakDetector

# Public API
__all__ = [
    'PeakModel',
    'PeakFitter',
    'PeakDetector'
    'PeakFitResult'
]