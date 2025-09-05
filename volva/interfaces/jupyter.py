"""
Volva XRD Analysis Tool - Jupyter Interface (SRM Enhanced)
Enhanced Jupyter notebook interface with SRM instrumental correction capabilities
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from typing import Dict, List, Optional, Union, Tuple
import logging

from ..core import (
    VolvaConfig, Constants, AnalysisError, ValidationError,
    LMFIT_AVAILABLE
)
from ..data import XRDDataProcessor
from ..peaks import PeakDetector
from ..analysis import (
    ScherrerAnalyzer, MonshiScherrerAnalyzer,
    WilliamsonHallAnalyzer, HalderWagnerAnalyzer
)
from ..srm import (
    SRMProcessor, InstrumentalCorrector,
    get_srm, list_available_srms
)

logger = logging.getLogger(__name__)


class JupyterVolva:
    """Enhanced Jupyter interface with SRM instrumental correction capabilities"""

    def __init__(self, wavelength_name: str = "Cu Ka", config_file: str = None):
        """
        Initialize Jupyter interface with SRM support

        Parameters:
        - wavelength_name: X-ray source name (from Constants.WAVELENGTHS)
        - config_file: Optional configuration file path
        """
        # Load configuration
        try:
            if config_file:
                self.config = VolvaConfig.load(config_file)
            else:
                self.config = VolvaConfig()
        except Exception as e:
            logger.warning(f"Using default config due to error: {e}")
            self.config = VolvaConfig()

        # Set up wavelength
        self.wavelength_name = wavelength_name
        if wavelength_name not in Constants.WAVELENGTHS:
            available = list(Constants.WAVELENGTHS.keys())
            raise ValueError(f"Unknown wavelength: {wavelength_name}. Available: {available}")

        self.wavelength = Constants.WAVELENGTHS[wavelength_name]['ka1']
        self.wavelength_info = Constants.WAVELENGTHS[wavelength_name]

        # Initialize components
        self.peak_detector = PeakDetector(self.config)
        self.scherrer_analyzer = ScherrerAnalyzer(self.wavelength, self.config)
        self.monshi_analyzer = MonshiScherrerAnalyzer(self.wavelength_info, self.config)
        self.wh_analyzer = WilliamsonHallAnalyzer(self.wavelength, self.config)
        self.hw_analyzer = HalderWagnerAnalyzer(self.wavelength, self.config)

        # SRM components
        self.srm_processor = SRMProcessor()
        self.instrumental_corrector = InstrumentalCorrector()

        # Data storage
        self.data = None
        self.srm_data = None
        self.peaks = None
        self.srm_result = None
        self.srm_calibration = None
        self.corrected_peaks = None
        self.last_results = {}

        # Display initialization info
        self._display_initialization_info()

    def _display_initialization_info(self):
        """Display initialization information with SRM features"""
        print("🔬 Volva XRD Analyzer v4.0 Enhanced (SRM) - Jupyter Interface")
        print("=" * 70)
        print(f"📡 Radiation: {self.wavelength_name} (λ = {self.wavelength:.5f} nm)")
        print(f"🧮 Advanced fitting: {'Available (lmfit)' if LMFIT_AVAILABLE else 'Basic only (install lmfit for advanced)'}")
        print(f"⚙️  Configuration: {len(self.config.instrumental_broadening_models)} instrument models available")
        print(f"📊 Analysis methods: Scherrer, Monshi-Scherrer, Williamson-Hall, Halder-Wagner")
        print(f"🔧 SRM Support: LaB₆, Silicon instrumental correction")
        print("=" * 70)
        print("💡 New SRM Workflow:")
        print("   1. load_srm_data() → 2. process_srm() → 3. create_calibration()")
        print("   4. load_data() → 5. detect_peaks() → 6. apply_corrections()")
        print("   7. run analysis methods with SRM-corrected peaks")
        print("=" * 70)

    def load_srm_data(self, file_path: str) -> bool:
        """
        Load SRM (Standard Reference Material) data

        Parameters:
        - file_path: Path to SRM data file

        Returns:
        - bool: Success status
        """
        try:
            two_theta, intensity, metadata = self.srm_processor.load_srm_data(file_path)

            self.srm_data = {
                'two_theta': two_theta,
                'intensity': intensity,
                'filename': os.path.basename(file_path),
                'full_path': file_path,
                'metadata': metadata
            }

            print(f"✅ SRM data loaded successfully")
            print(f"   📊 File: {os.path.basename(file_path)}")
            print(f"   📈 Range: {two_theta.min():.1f}° - {two_theta.max():.1f}°")
            print(f"   📊 Points: {len(two_theta):,}")
            print(f"   💡 Next: Use process_srm() to identify and process SRM peaks")

            return True

        except Exception as e:
            print(f"❌ SRM data loading failed: {str(e)}")
            return False

    def identify_srm_material(self, candidate_materials: Optional[List[str]] = None) -> Optional[str]:
        """
        Automatically identify SRM material from loaded data

        Parameters:
        - candidate_materials: List of materials to consider (all if None)

        Returns:
        - Detected SRM material name or None
        """
        if self.srm_data is None:
            print("❌ Please load SRM data first using load_srm_data()")
            return None

        try:
            print("🔍 Auto-detecting SRM material...")

            detected_material = self.srm_processor.identify_srm_material(
                self.srm_data['two_theta'],
                self.srm_data['intensity'],
                self.wavelength,
                candidate_materials
            )

            print(f"✅ SRM material detected: {detected_material}")
            print(f"   💡 Next: Use process_srm('{detected_material}') to process peaks")

            return detected_material

        except Exception as e:
            print(f"❌ SRM material identification failed: {str(e)}")
            return None

    def process_srm(self, srm_material: str = None, position_tolerance: float = 0.15) -> bool:
        """
        Process SRM data to identify and match theoretical peaks

        Parameters:
        - srm_material: SRM material name (auto-detect if None)
        - position_tolerance: Peak matching tolerance in degrees

        Returns:
        - bool: Success status
        """
        if self.srm_data is None:
            print("❌ Please load SRM data first using load_srm_data()")
            return False

        try:
            # Auto-detect material if not specified
            if srm_material is None:
                srm_material = self.identify_srm_material()
                if srm_material is None:
                    return False

            print(f"🔧 Processing SRM peaks for {srm_material}...")

            # Process SRM data
            self.srm_result = self.srm_processor.process_srm_data(
                self.srm_data['two_theta'],
                self.srm_data['intensity'],
                srm_material,
                self.wavelength,
                position_tolerance=position_tolerance
            )

            # Display results
            matched = len(self.srm_result.matched_peaks)
            missing = len(self.srm_result.missing_expected)
            unmatched = len(self.srm_result.unmatched_peaks)
            quality = self.srm_result.quality_score

            print(f"✅ SRM processing complete:")
            print(f"   📊 Material: {self.srm_result.srm_name}")
            print(f"   ✅ Matched peaks: {matched}")
            print(f"   ❌ Missing expected: {missing}")
            print(f"   ⚠️  Unmatched peaks: {unmatched}")
            print(f"   📈 Quality score: {quality:.3f}")

            if quality < 0.7:
                print(f"   ⚠️  Warning: Low quality score - check SRM data")

            print(f"   💡 Next: Use create_calibration() to fit instrumental profile")

            return True

        except Exception as e:
            print(f"❌ SRM processing failed: {str(e)}")
            return False

    def create_calibration(self, function_name: str = None, auto_select: bool = True) -> bool:
        """
        Create instrumental calibration from processed SRM data

        Parameters:
        - function_name: Instrumental function ('Caglioti', 'Polynomial', 'Pearson-VII')
        - auto_select: Automatically select best function

        Returns:
        - bool: Success status
        """
        if self.srm_result is None:
            print("❌ Please process SRM data first using process_srm()")
            return False

        try:
            print("🔧 Creating instrumental calibration...")

            # Create instrument configuration
            instrument_config = {
                'voltage': 40,  # kV
                'current': 40,  # mA
                'detector': 'Standard',
                'wavelength': self.wavelength,
                'radiation': self.wavelength_name
            }

            # Create calibration
            self.srm_calibration = self.instrumental_corrector.create_calibration(
                self.srm_result,
                instrument_config,
                function_name=function_name,
                auto_select_function=auto_select
            )

            # Display results
            func_name = self.srm_calibration.fit_result.function_name
            r_squared = self.srm_calibration.fit_result.r_squared
            validity_range = self.srm_calibration.validity_range

            print(f"✅ Instrumental calibration created:")
            print(f"   📊 Function: {func_name}")
            print(f"   📈 R²: {r_squared:.4f}")
            print(f"   📐 Valid range: {validity_range[0]:.1f}° - {validity_range[1]:.1f}°")
            print(f"   📋 Parameters: {len(self.srm_calibration.fit_result.parameters)}")

            # Display parameters
            for param, value in self.srm_calibration.fit_result.parameters.items():
                print(f"      {param}: {value:.6f}")

            if r_squared < 0.9:
                print(f"   ⚠️  Warning: R² < 0.9 - consider different function or check data")

            print(f"   💡 Next: Load sample data and apply corrections with apply_corrections()")

            return True

        except Exception as e:
            print(f"❌ Calibration creation failed: {str(e)}")
            return False

    def save_calibration(self, filename: str) -> bool:
        """Save SRM calibration to file"""
        if self.srm_calibration is None:
            print("❌ No calibration to save")
            return False

        try:
            self.instrumental_corrector.save_calibration(self.srm_calibration, filename)
            print(f"✅ Calibration saved to {filename}")
            return True
        except Exception as e:
            print(f"❌ Failed to save calibration: {str(e)}")
            return False

    def load_calibration(self, filename: str) -> bool:
        """Load SRM calibration from file"""
        try:
            self.srm_calibration = self.instrumental_corrector.load_calibration(filename)
            func_name = self.srm_calibration.fit_result.function_name
            r_squared = self.srm_calibration.fit_result.r_squared
            print(f"✅ Calibration loaded: {func_name} (R² = {r_squared:.4f})")
            return True
        except Exception as e:
            print(f"❌ Failed to load calibration: {str(e)}")
            return False

    def apply_corrections(self, method: str = 'deconvolution') -> bool:
        """
        Apply SRM corrections to detected sample peaks

        Parameters:
        - method: Correction method ('deconvolution', 'subtraction', 'both')

        Returns:
        - bool: Success status
        """
        if self.peaks is None:
            print("❌ Please detect sample peaks first using detect_peaks()")
            return False

        if self.srm_calibration is None:
            print("❌ Please create SRM calibration first using create_calibration()")
            return False

        try:
            print(f"🔧 Applying SRM corrections using {method} method...")

            # Convert peaks to PeakFitResult format
            from ..peaks.models import PeakFitResult
            peak_results = []
            for peak in self.peaks:
                peak_result = PeakFitResult(
                    position=peak['2theta'],
                    intensity=peak['intensity'],
                    fwhm=peak['FWHM_obs'],
                    area=peak.get('area', peak['intensity'] * peak['FWHM_obs']),
                    background=peak.get('background', 0),
                    r_squared=peak.get('r_squared', 0.85)
                )
                peak_results.append(peak_result)

            # Apply corrections
            self.corrected_peaks = self.instrumental_corrector.apply_correction(
                peak_results,
                self.srm_calibration,
                method=method
            )

            # Assess correction quality
            correction_quality = self.instrumental_corrector.estimate_correction_quality(
                self.corrected_peaks
            )

            # Display results
            print(f"✅ SRM corrections applied:")
            print(f"   📊 Method: {method}")
            print(f"   📈 Peaks corrected: {len(self.corrected_peaks)}")
            print(f"   📊 Quality: {correction_quality['quality']}")

            metrics = correction_quality.get('metrics', {})
            if 'deconvolution_rate' in metrics:
                print(f"   ✅ Valid deconvolutions: {metrics['deconvolution_rate']:.1%}")
            if 'mean_correction_factor' in metrics:
                print(f"   📉 Mean correction factor: {metrics['mean_correction_factor']:.3f}")

            # Show correction summary
            print(f"\n   📋 Peak Correction Summary:")
            for i, (orig, corr) in enumerate(zip(self.peaks[:5], self.corrected_peaks[:5])):
                ratio = corr.corrected_fwhm / corr.observed_fwhm
                status = "✅" if corr.deconvolution_valid else "⚠️"
                print(f"      Peak {i+1}: {corr.observed_fwhm:.3f}° → {corr.corrected_fwhm:.3f}° ({ratio:.2f}) {status}")

            if len(self.corrected_peaks) > 5:
                print(f"      ... and {len(self.corrected_peaks)-5} more peaks")

            print(f"   💡 Next: Run analysis methods - they will automatically use corrected peaks")

            return True

        except Exception as e:
            print(f"❌ SRM correction failed: {str(e)}")
            return False

    # Enhanced load_data with SRM status
    def load_data(self, file_path: str) -> bool:
        """
        Load XRD data from file with enhanced error handling

        Parameters:
        - file_path: Path to XRD data file

        Returns:
        - bool: Success status
        """
        try:
            two_theta, intensity, message = XRDDataProcessor.load_xrd_file(file_path)

            self.data = {
                'two_theta': two_theta,
                'intensity': intensity,
                'filename': os.path.basename(file_path),
                'full_path': file_path
            }

            print(f"✅ {message}")
            print(f"   📊 Range: {two_theta.min():.1f}° - {two_theta.max():.1f}°")
            print(f"   📈 Intensity: {intensity.min():.0f} - {intensity.max():.0f} counts")

            # Validate data quality
            if XRDDataProcessor.validate_xrd_data(two_theta, intensity):
                print("   ✅ Data validation passed")
            else:
                print("   ⚠️  Data validation warnings - check data quality")

            # Show SRM status
            if self.srm_calibration:
                print(f"   🔧 SRM calibration available: {self.srm_calibration.srm_name}")
            else:
                print(f"   💡 Tip: Load SRM data first for best accuracy")

            self.plot_data()
            return True

        except Exception as e:
            print(f"❌ Data loading failed: {str(e)}")
            return False

    # Enhanced analysis methods with SRM support
    def run_scherrer(self, k_constant: Union[str, float] = None,
                     filter_quality: bool = True, use_srm: bool = None) -> Optional[Dict]:
        """
        Enhanced Scherrer analysis with SRM correction support

        Parameters:
        - k_constant: Scherrer shape factor
        - filter_quality: Filter peaks by quality before analysis
        - use_srm: Force SRM usage (auto-detect if None)

        Returns:
        - Analysis results dictionary
        """
        if self.peaks is None:
            print("❌ Please detect peaks first using detect_peaks()")
            return None

        try:
            peaks_to_analyze = self.peaks
            if filter_quality:
                peaks_to_analyze = self.peak_detector.filter_peaks_by_quality(
                    self.peaks, self.config.min_r_squared)
                if len(peaks_to_analyze) < len(self.peaks):
                    print(f"🔍 Quality filtering: {len(peaks_to_analyze)}/{len(self.peaks)} peaks retained")

            # Determine SRM usage
            using_srm = self.corrected_peaks is not None
            if use_srm is not None:
                using_srm = use_srm and self.corrected_peaks is not None

            print(f"🧮 Running Scherrer analysis{'with SRM corrections' if using_srm else ''}...")

            if using_srm:
                result = self.scherrer_analyzer.analyze_with_srm(
                    peaks_to_analyze, self.corrected_peaks, self.srm_calibration, k_constant)
            else:
                result = self.scherrer_analyzer.analyze(peaks_to_analyze, k_constant)

            if result['success']:
                self._display_scherrer_results(result)
                self.last_results['scherrer'] = result
                return result
            else:
                print(f"❌ Scherrer analysis failed: {result.get('message', 'Unknown error')}")
                return None

        except Exception as e:
            print(f"❌ Scherrer analysis failed: {str(e)}")
            return None

    def _display_scherrer_results(self, result: Dict):
        """Display Scherrer analysis results with SRM info"""
        print(f"\n📊 {result['method']} Results:")
        print(f"   💎 Average Size: {result['average_size']:.1f} ± {result['std_size']:.1f} nm")
        print(f"   📊 Median Size: {result['median_size']:.1f} nm")
        print(f"   🎯 Weighted Average: {result.get('weighted_average_size', result['average_size']):.1f} nm")
        print(f"   📢 Peaks Used: {result['num_peaks']}")
        print(f"   ⚙️  K Constant: {result['k_constant']:.3f} ({result.get('k_description', '')})")

        # SRM information
        srm_info = result.get('srm_info', {})
        if srm_info.get('srm_corrected', False):
            print(f"   🔧 SRM Corrected: YES ({srm_info.get('srm_material', 'Unknown')})")
            quality = srm_info.get('correction_quality', {}).get('quality', 'unknown')
            print(f"   📈 Correction Quality: {quality}")
        else:
            print(f"   🔧 SRM Corrected: NO")

        if result.get('failed_peaks'):
            print(f"   ⚠️  Failed Peaks: {len(result['failed_peaks'])}")

        # Display quality assessment
        assessment = result.get('assessment', {})
        quality = assessment.get('overall_quality', 'unknown')
        confidence = assessment.get('confidence_level', 'unknown')
        print(f"   🎭 Quality: {quality} (confidence: {confidence})")

        # SRM recommendations
        srm_recommendations = srm_info.get('srm_recommendations', [])
        srm_warnings = srm_info.get('srm_warnings', [])

        if srm_recommendations:
            print(f"   💡 SRM Recommendations:")
            for rec in srm_recommendations[:2]:  # Show first 2
                print(f"      • {rec}")

        if srm_warnings:
            print(f"   ⚠️  SRM Warnings:")
            for warning in srm_warnings[:2]:  # Show first 2
                print(f"      • {warning}")

    # Similar updates for other analysis methods
    def run_monshi_scherrer(self, k_constant: Union[str, float] = None,
                            doublet_method: str = 'correction',
                            filter_quality: bool = True, use_srm: bool = None) -> Optional[Dict]:
        """Enhanced Monshi-Scherrer analysis with SRM support"""
        if self.peaks is None:
            print("❌ Please detect peaks first using detect_peaks()")
            return None

        try:
            peaks_to_analyze = self.peaks
            if filter_quality:
                peaks_to_analyze = self.peak_detector.filter_peaks_by_quality(
                    self.peaks, self.config.min_r_squared)

            using_srm = self.corrected_peaks is not None
            if use_srm is not None:
                using_srm = use_srm and self.corrected_peaks is not None

            print(f"🧮 Running Monshi-Scherrer analysis{'with SRM corrections' if using_srm else ''} (doublet: {doublet_method})...")

            if using_srm:
                result = self.monshi_analyzer.analyze_with_srm(
                    peaks_to_analyze, self.corrected_peaks, self.srm_calibration,
                    k_constant, doublet_method)
            else:
                result = self.monshi_analyzer.analyze(peaks_to_analyze, k_constant, doublet_method)

            if result['success']:
                self._display_monshi_results(result)
                self.last_results['monshi'] = result
                return result
            else:
                print(f"❌ Monshi-Scherrer analysis failed")
                return None

        except Exception as e:
            print(f"❌ Monshi-Scherrer analysis failed: {str(e)}")
            return None

    def _display_monshi_results(self, result: Dict):
        """Display Monshi-Scherrer results with SRM info"""
        print(f"\n📊 {result['method']} Results:")
        print(f"   💎 Average Size: {result['average_size']:.1f} ± {result['std_size']:.1f} nm")
        print(f"   🎯 Weighted Average: {result.get('weighted_average_size', result['average_size']):.1f} nm")
        print(f"   📢 Peaks Used: {result['num_peaks']}")
        print(f"   📡 Kα1: {result['ka1_wavelength']:.5f} nm, Kα2: {result['ka2_wavelength']:.5f} nm")
        print(f"   🔧 Alpha Parameter: {result['alpha_parameter']:.3e}")
        print(f"   📝 Doublet Method: {result['doublet_method']}")

        # SRM information
        srm_info = result.get('srm_info', {})
        if srm_info.get('srm_corrected', False):
            print(f"   🔧 SRM Corrected: YES ({srm_info.get('srm_material', 'Unknown')})")
        else:
            print(f"   🔧 SRM Corrected: NO")

    def run_williamson_hall(self, k_constant: float = None,
                            use_robust_fitting: bool = True,
                            filter_quality: bool = True, use_srm: bool = None) -> Optional[Dict]:
        """Enhanced Williamson-Hall analysis with SRM support"""
        if self.peaks is None:
            print("❌ Please detect peaks first using detect_peaks()")
            return None

        peaks_to_analyze = self.peaks
        if filter_quality:
            peaks_to_analyze = self.peak_detector.filter_peaks_by_quality(
                self.peaks, self.config.min_r_squared)

        if len(peaks_to_analyze) < 3:
            print(f"❌ Need at least 3 peaks for W-H analysis. Found {len(peaks_to_analyze)}")
            return None

        try:
            using_srm = self.corrected_peaks is not None
            if use_srm is not None:
                using_srm = use_srm and self.corrected_peaks is not None

            print(f"🧮 Running Williamson-Hall analysis{'with SRM corrections' if using_srm else ''}...")

            if using_srm:
                result = self.wh_analyzer.analyze_with_srm(
                    peaks_to_analyze, self.corrected_peaks, self.srm_calibration,
                    k_constant, use_robust_fitting)
            else:
                result = self.wh_analyzer.analyze(peaks_to_analyze, k_constant, use_robust_fitting)

            if result['success']:
                self._display_williamson_hall_results(result)
                self.plot_williamson_hall(result)
                self.last_results['williamson_hall'] = result
                return result
            else:
                print(f"❌ W-H analysis failed")
                return None

        except Exception as e:
            print(f"❌ W-H analysis failed: {str(e)}")
            return None

    def _display_williamson_hall_results(self, result: Dict):
        """Display Williamson-Hall results with SRM info"""
        print(f"\n📊 {result['method']} Results:")
        print(f"   💎 Crystallite Size: {result['crystallite_size_nm']:.1f} ± {result.get('crystallite_size_uncertainty_nm', 0):.1f} nm")
        print(f"   🔧 Microstrain: {result['microstrain_percent']:.5f} ± {result.get('microstrain_uncertainty_percent', 0):.5f} %")
        print(f"   📊 Linear Fit R²: {result['r_squared']:.4f}")
        print(f"   📢 Peaks Used: {result['num_peaks']}")

        # SRM information
        srm_info = result.get('srm_info', {})
        if srm_info.get('srm_corrected', False):
            print(f"   🔧 SRM Corrected: YES ({srm_info.get('srm_material', 'Unknown')})")
            quality = srm_info.get('correction_quality', {}).get('quality', 'unknown')
            print(f"   📈 Correction Quality: {quality}")
        else:
            print(f"   🔧 SRM Corrected: NO")

        # Display isotropy assessment
        isotropy = result.get('isotropy_assessment', {})
        if isotropy:
            is_isotropic = isotropy.get('likely_isotropic', False)
            confidence = isotropy.get('confidence', 'unknown')
            print(f"   🎯 Isotropy: {'✅ Likely' if is_isotropic else '❌ Questionable'} (confidence: {confidence})")

        # Overall assessment
        assessment = result.get('assessment', {})
        quality = assessment.get('overall_quality', 'unknown')
        validity = assessment.get('scientific_validity', 'unknown')
        print(f"   🎭 Overall Quality: {quality}")
        print(f"   🔬 Scientific Validity: {validity}")

        def run_halder_wagner(self, a_parameter: float = None, k_constant: float = None,
                              use_robust_fitting: bool = True, compare_with_wh: bool = True,
                              filter_quality: bool = True, use_srm: bool = None) -> Optional[Dict]:
            """Enhanced Halder-Wagner analysis with SRM support"""
            if self.peaks is None:
                print("Please detect peaks first using detect_peaks()")
                return None

            peaks_to_analyze = self.peaks
            if filter_quality:
                peaks_to_analyze = self.peak_detector.filter_peaks_by_quality(
                    self.peaks, self.config.min_r_squared)

            if len(peaks_to_analyze) < 3:
                print(f"Need at least 3 peaks for H-W analysis. Found {len(peaks_to_analyze)}")
                return None

            try:
                using_srm = self.corrected_peaks is not None
                if use_srm is not None:
                    using_srm = use_srm and self.corrected_peaks is not None

                print(f"Running Halder-Wagner analysis{' with SRM corrections' if using_srm else ''}...")

                if using_srm:
                    result = self.hw_analyzer.analyze_with_srm(
                        peaks_to_analyze, self.corrected_peaks, self.srm_calibration,
                        a_parameter, k_constant, use_robust_fitting, compare_with_wh)
                else:
                    result = self.hw_analyzer.analyze(
                        peaks_to_analyze, a_parameter, k_constant,
                        use_robust_fitting, compare_with_wh)

                if result['success']:
                    self._display_halder_wagner_results(result)
                    self.plot_halder_wagner(result)
                    self.last_results['halder_wagner'] = result
                    return result
                else:
                    print("H-W analysis failed")
                    return None

            except Exception as e:
                print(f"H-W analysis failed: {str(e)}")
                return None

        def _display_halder_wagner_results(self, result: Dict):
            """Display Halder-Wagner results with SRM info"""
            print(f"\n{result['method']} Results:")
            print(
                f"   Crystallite Size: {result['crystallite_size_nm']:.1f} ± {result.get('crystallite_size_uncertainty_nm', 0):.1f} nm")
            print(
                f"   Microstrain: {result['microstrain_percent']:.5f} ± {result.get('microstrain_uncertainty_percent', 0):.5f} %")
            print(f"   Fit R2: {result['r_squared']:.4f}")
            print(f"   Peaks Used: {result['num_peaks']}")
            print(f"   'a' Parameter: {result['a_parameter']:.1f}")

            # Display model information
            quality_metrics = result.get('quality_metrics', {})
            model_type = quality_metrics.get('model_type', 'unknown')
            print(f"   Model Type: {model_type}")

            # SRM information
            srm_info = result.get('srm_info', {})
            if srm_info.get('srm_corrected', False):
                print(f"   SRM Corrected: YES ({srm_info.get('srm_material', 'Unknown')})")
                quality = srm_info.get('correction_quality', {}).get('quality', 'unknown')
                print(f"   Correction Quality: {quality}")
            else:
                print(f"   SRM Corrected: NO")

            # Display anisotropy assessment
            anisotropy = result.get('anisotropy_assessment', {})
            if anisotropy:
                is_anisotropic = anisotropy.get('likely_anisotropic', False)
                confidence = anisotropy.get('confidence', 'unknown')
                print(
                    f"   Anisotropy: {'Detected' if is_anisotropic else 'Limited evidence'} (confidence: {confidence})")

            # W-H comparison
            wh_comparison = result.get('wh_comparison', {})
            if wh_comparison.get('comparison_performed', False):
                wh_r2 = wh_comparison.get('wh_r_squared', 0)
                improvement = result['r_squared'] - wh_r2
                print(f"   vs W-H: R2 improvement = {improvement:+.3f}")

        def run_all_methods(self, verbose: bool = True, compare_results: bool = True,
                            use_srm: bool = None) -> Dict[str, Dict]:
            """
            Run all analysis methods with comprehensive comparison including SRM support

            Parameters:
            - verbose: Display detailed output for each method
            - compare_results: Generate comparison analysis
            - use_srm: Force SRM usage (auto-detect if None)

            Returns:
            - Dictionary of all results
            """
            if self.peaks is None:
                print("Please detect peaks first using detect_peaks()")
                return {}

            print("Running comprehensive analysis with all methods...")
            print("=" * 60)

            # Show SRM status
            using_srm = self.corrected_peaks is not None
            if use_srm is not None:
                using_srm = use_srm and self.corrected_peaks is not None

            if using_srm:
                print(f"SRM Corrections: ACTIVE ({self.srm_calibration.srm_name})")
            else:
                print("SRM Corrections: NOT APPLIED (traditional method)")
            print("=" * 60)

            results = {}

            # Run size-only methods
            print("\nSIZE-ONLY METHODS")
            print("-" * 30)

            scherrer_result = self.run_scherrer(use_srm=use_srm)
            if scherrer_result:
                results['scherrer'] = scherrer_result

            monshi_result = self.run_monshi_scherrer(use_srm=use_srm)
            if monshi_result:
                results['monshi'] = monshi_result

            # Run size + strain methods
            if len(self.peaks) >= 3:
                print("\nSIZE + STRAIN METHODS")
                print("-" * 30)

                wh_result = self.run_williamson_hall(use_srm=use_srm)
                if wh_result:
                    results['williamson_hall'] = wh_result

                hw_result = self.run_halder_wagner(use_srm=use_srm)
                if hw_result:
                    results['halder_wagner'] = hw_result
            else:
                print("\nInsufficient peaks for strain analysis methods (need ≥3)")

            # Generate comparison if requested
            if compare_results and len(results) >= 2:
                print("\nMETHOD COMPARISON")
                print("-" * 30)
                self._display_method_comparison(results)

                if len(results) >= 2:
                    self.plot_method_comparison(results)

            self.last_results = results
            print(f"\nAnalysis complete! {len(results)} methods successful")

            # SRM summary
            if using_srm:
                print("Results include SRM instrumental broadening correction")
            else:
                print("Consider using SRM correction for higher accuracy")

            return results

        def _display_method_comparison(self, results: Dict[str, Dict]):
            """Display comparison between analysis methods with SRM awareness"""
            # Size comparison
            size_methods = {}
            strain_methods = {}

            for method, result in results.items():
                if 'average_size' in result:  # Size-only methods
                    size_methods[method] = {
                        'name': result['method'],
                        'size': result['average_size'],
                        'std': result['std_size'],
                        'peaks': result['num_peaks'],
                        'srm': result.get('srm_info', {}).get('srm_corrected', False)
                    }
                elif 'crystallite_size_nm' in result:  # Size + strain methods
                    strain_methods[method] = {
                        'name': result['method'],
                        'size': result['crystallite_size_nm'],
                        'strain': result['microstrain_percent'],
                        'r2': result['r_squared'],
                        'peaks': result['num_peaks'],
                        'srm': result.get('srm_info', {}).get('srm_corrected', False)
                    }

            if size_methods:
                print("SIZE METHODS COMPARISON:")
                for method, data in size_methods.items():
                    srm_status = "SRM" if data['srm'] else "STD"
                    print(
                        f"   {data['name']:20}: {data['size']:6.1f} ± {data['std']:4.1f} nm (n={data['peaks']}) [{srm_status}]")

                if len(size_methods) > 1:
                    sizes = [data['size'] for data in size_methods.values()]
                    cv = np.std(sizes) / np.mean(sizes) * 100
                    print(f"   {'Coefficient of Variation':<20}: {cv:6.1f}%")

            if strain_methods:
                print("\nSIZE + STRAIN METHODS:")
                print(f"   {'Method':<20} {'Size (nm)':<10} {'Strain (%)':<12} {'R2':<8} {'Peaks':<6} {'SRM':<4}")
                print("   " + "-" * 66)
                for method, data in strain_methods.items():
                    srm_status = "YES" if data['srm'] else "NO"
                    print(f"   {data['name']:<20} {data['size']:<10.1f} "
                          f"{data['strain']:<12.5f} {data['r2']:<8.4f} {data['peaks']:<6} {srm_status}")

            # Method recommendations with SRM awareness
            print("\nRECOMMENDATIONS:")
            self._generate_method_recommendations(results)

        def _generate_method_recommendations(self, results: Dict[str, Dict]):
            """Generate recommendations based on analysis results with SRM awareness"""
            recommendations = []

            # Check SRM usage consistency
            srm_usage = [r.get('srm_info', {}).get('srm_corrected', False) for r in results.values()]
            if any(srm_usage) and not all(srm_usage):
                recommendations.append("Mixed SRM usage detected - use consistently for best comparison")
            elif all(srm_usage):
                recommendations.append("Excellent: All methods use SRM correction")
            elif not any(srm_usage) and self.srm_calibration:
                recommendations.append("SRM calibration available but not used - apply corrections for best accuracy")

            # Check for method agreement
            size_methods = [r for r in results.values() if 'average_size' in r]
            if len(size_methods) >= 2:
                sizes = [r['average_size'] for r in size_methods]
                cv = np.std(sizes) / np.mean(sizes) * 100
                if cv < 10:
                    recommendations.append("Size methods show excellent agreement")
                elif cv < 20:
                    recommendations.append("Size methods show moderate agreement")
                else:
                    recommendations.append("Size methods show poor agreement - check sample homogeneity")

            # Check strain analysis quality
            strain_methods = [r for r in results.values() if 'microstrain_percent' in r]
            for result in strain_methods:
                method_name = result['method']
                r2 = result['r_squared']

                if method_name == 'Williamson-Hall':
                    isotropy = result.get('isotropy_assessment', {})
                    if not isotropy.get('likely_isotropic', True):
                        recommendations.append("W-H isotropy assumption questionable - consider H-W method")
                elif method_name == 'Halder-Wagner':
                    anisotropy = result.get('anisotropy_assessment', {})
                    if anisotropy.get('likely_anisotropic', False):
                        recommendations.append("H-W method appropriate for detected anisotropic strain")

            # Display recommendations
            for rec in recommendations:
                print(f"   • {rec}")

        # Enhanced plotting methods
        def plot_data(self):
            """Plot XRD data with SRM status"""
            if self.data is None:
                print("No data to plot")
                return

            fig, ax = plt.subplots(figsize=(12, 6))
            ax.plot(self.data['two_theta'], self.data['intensity'],
                    'b-', linewidth=1.5, alpha=0.8, label='Sample Data')

            ax.set_xlabel('2θ (degrees)', fontsize=12)
            ax.set_ylabel('Intensity (counts)', fontsize=12)

            title = f'XRD Pattern: {self.data["filename"]} ({self.wavelength_name})'
            if self.srm_calibration:
                title += f' [SRM: {self.srm_calibration.srm_name}]'
            ax.set_title(title, fontsize=13, fontweight='bold')

            ax.grid(True, alpha=0.3)
            ax.legend()

            plt.tight_layout()
            plt.show()

        def plot_data_with_peaks(self):
            """Plot XRD data with detected peaks and SRM corrections"""
            if self.data is None or self.peaks is None:
                print("No data or peaks to plot")
                return

            fig, ax = plt.subplots(figsize=(12, 6))
            ax.plot(self.data['two_theta'], self.data['intensity'],
                    'b-', linewidth=1.5, alpha=0.8, label='Sample Data')

            # Plot peaks with quality-based coloring
            peak_positions = [p['2theta'] for p in self.peaks]
            peak_intensities = [p['intensity'] for p in self.peaks]
            peak_qualities = [p.get('r_squared', 0.5) for p in self.peaks]

            colors = ['red' if q > 0.8 else 'orange' if q > 0.6 else 'yellow' for q in peak_qualities]

            ax.scatter(peak_positions, peak_intensities, c=colors, s=100,
                       marker='^', zorder=5, alpha=0.9, edgecolors='black',
                       linewidth=1, label='Detected Peaks')

            # Peak labels with SRM correction info
            for i, (pos, intensity, fwhm) in enumerate(zip(
                    peak_positions, peak_intensities, [p['FWHM_obs'] for p in self.peaks])):

                if self.corrected_peaks and i < len(self.corrected_peaks):
                    corrected_fwhm = self.corrected_peaks[i].corrected_fwhm
                    label = f'P{i + 1}\nObs: {fwhm:.3f}°\nCorr: {corrected_fwhm:.3f}°'
                else:
                    label = f'P{i + 1}\n{fwhm:.3f}°'

                ax.annotate(label, (pos, intensity),
                            xytext=(5, 15), textcoords='offset points',
                            fontsize=8, ha='center',
                            bbox=dict(boxstyle="round,pad=0.2", facecolor=colors[i], alpha=0.3))

            ax.set_xlabel('2θ (degrees)', fontsize=12)
            ax.set_ylabel('Intensity (counts)', fontsize=12)

            title = f'Peak Detection: {self.data["filename"]} ({len(self.peaks)} peaks)'
            if self.corrected_peaks:
                title += ' [SRM Corrected]'
            ax.set_title(title, fontsize=13, fontweight='bold')

            ax.grid(True, alpha=0.3)
            ax.legend()

            plt.tight_layout()
            plt.show()

        # Keep existing plotting methods from original (plot_williamson_hall, plot_halder_wagner, etc.)
        # but add SRM info to titles and annotations

        def show_srm_status(self):
            """Display current SRM workflow status"""
            print("\nSRM WORKFLOW STATUS")
            print("=" * 50)

            # SRM Data
            if self.srm_data:
                print(f"1. SRM Data: LOADED ({self.srm_data['filename']})")
            else:
                print("1. SRM Data: NOT LOADED")

            # SRM Processing
            if self.srm_result:
                matched = len(self.srm_result.matched_peaks)
                quality = self.srm_result.quality_score
                print(
                    f"2. SRM Processing: COMPLETE ({self.srm_result.srm_name}, {matched} peaks, quality: {quality:.3f})")
            else:
                print("2. SRM Processing: NOT DONE")

            # Calibration
            if self.srm_calibration:
                func_name = self.srm_calibration.fit_result.function_name
                r2 = self.srm_calibration.fit_result.r_squared
                print(f"3. Calibration: CREATED ({func_name}, R²: {r2:.4f})")
            else:
                print("3. Calibration: NOT CREATED")

            # Sample Data
            if self.data:
                print(f"4. Sample Data: LOADED ({self.data['filename']})")
            else:
                print("4. Sample Data: NOT LOADED")

            # Peak Detection
            if self.peaks:
                print(f"5. Peak Detection: COMPLETE ({len(self.peaks)} peaks)")
            else:
                print("5. Peak Detection: NOT DONE")

            # Corrections
            if self.corrected_peaks:
                quality = "Good" if len([p for p in self.corrected_peaks if p.deconvolution_valid]) > len(
                    self.corrected_peaks) * 0.8 else "Fair"
                print(f"6. SRM Corrections: APPLIED ({len(self.corrected_peaks)} peaks, {quality})")
            else:
                print("6. SRM Corrections: NOT APPLIED")

            print("=" * 50)

            # Next steps
            if not self.srm_data:
                print("Next: load_srm_data('path/to/srm_file.dat')")
            elif not self.srm_result:
                print("Next: process_srm() or identify_srm_material()")
            elif not self.srm_calibration:
                print("Next: create_calibration()")
            elif not self.data:
                print("Next: load_data('path/to/sample_file.dat')")
            elif not self.peaks:
                print("Next: detect_peaks()")
            elif not self.corrected_peaks:
                print("Next: apply_corrections()")
            else:
                print("Ready: Run analysis methods (they will use SRM corrections)")

        def show_help(self):
            """Display comprehensive help with SRM workflow"""
            help_text = """
    VOLVA SRM-ENHANCED JUPYTER INTERFACE HELP
    ==========================================

    SRM WORKFLOW (Recommended):
    1. load_srm_data('srm_file.dat') - Load LaB6 or Silicon standard
    2. identify_srm_material() - Auto-detect material (optional)
    3. process_srm('LaB6') - Process and match SRM peaks
    4. create_calibration() - Fit instrumental broadening function
    5. load_data('sample_file.dat') - Load your sample data
    6. detect_peaks() - Find sample peaks
    7. apply_corrections() - Apply SRM corrections to peaks
    8. Run analysis methods - Automatically use corrected peaks

    TRADITIONAL WORKFLOW (Basic):
    1. load_data('sample_file.dat')
    2. detect_peaks()
    3. Run analysis methods - Uses traditional instrumental correction

    KEY METHODS:
    • load_srm_data(file_path) - Load SRM data
    • process_srm(material) - Process SRM peaks
    • create_calibration() - Create instrumental profile
    • apply_corrections() - Apply SRM corrections
    • load_data(file_path) - Load sample data
    • detect_peaks() - Find peaks with quality assessment
    • run_scherrer() - Basic crystallite size
    • run_monshi_scherrer() - Enhanced with Kα2 correction
    • run_williamson_hall() - Size + strain with isotropy check
    • run_halder_wagner() - Anisotropic strain analysis
    • run_all_methods() - Comprehensive analysis with comparison

    SRM UTILITIES:
    • show_srm_status() - Check workflow progress
    • save_calibration(filename) - Save calibration for reuse
    • load_calibration(filename) - Load existing calibration
    • list_available_srms() - Show supported SRM materials

    BENEFITS OF SRM CORRECTION:
    • Accurate crystallite sizes (removes instrumental broadening)
    • Reliable microstrain values
    • Angle-dependent correction (better than fixed values)
    • Publication-quality results

    For detailed parameter info, use help(method_name)
            """
            print(help_text)

        # Keep all other existing methods from the original jupyter.py
        # (export_results, save_config, load_config, get_peak_summary, etc.)
        # but ensure they're compatible with SRM data structures

        def export_results(self, filename: str, format_type: str = 'json') -> bool:
            """Export analysis results to file with SRM information"""
            if not self.last_results:
                print("No analysis results to export")
                return False

            try:
                if format_type.lower() == 'json':
                    export_data = {
                        'metadata': {
                            'software': 'Volva v4.0 Enhanced (SRM) - Jupyter Interface',
                            'analysis_date': datetime.now().isoformat(),
                            'sample_filename': self.data.get('filename') if self.data else None,
                            'wavelength_name': self.wavelength_name,
                            'wavelength_nm': self.wavelength,
                            'num_peaks': len(self.peaks) if self.peaks else 0,
                            'srm_corrected': self.corrected_peaks is not None,
                            'srm_material': self.srm_calibration.srm_name if self.srm_calibration else None
                        },
                        'srm_data': {
                            'calibration': {
                                'function': self.srm_calibration.fit_result.function_name if self.srm_calibration else None,
                                'r_squared': self.srm_calibration.fit_result.r_squared if self.srm_calibration else None,
                                'parameters': self.srm_calibration.fit_result.parameters if self.srm_calibration else None
                            } if self.srm_calibration else None,
                            'correction_applied': self.corrected_peaks is not None
                        },
                        'configuration': self.config.to_dict(),
                        'peaks': self.peaks if self.peaks else [],
                        'results': self.last_results
                    }

                    # Convert numpy types for JSON serialization
                    def convert_numpy(obj):
                        if isinstance(obj, np.ndarray):
                            return obj.tolist()
                        elif isinstance(obj, (np.floating, np.complexfloating)):
                            return float(obj)
                        elif isinstance(obj, (np.integer)):
                            return int(obj)
                        return obj

                    def clean_for_json(obj):
                        if isinstance(obj, dict):
                            return {k: clean_for_json(v) for k, v in obj.items()}
                        elif isinstance(obj, list):
                            return [clean_for_json(v) for v in obj]
                        else:
                            return convert_numpy(obj)

                    clean_export_data = clean_for_json(export_data)

                    with open(filename, 'w') as f:
                        json.dump(clean_export_data, f, indent=2)

                else:
                    raise ValueError(f"Unsupported format: {format_type}")

                print(f"Results exported to {filename}")
                return True

            except Exception as e:
                print(f"Export failed: {str(e)}")
                return False