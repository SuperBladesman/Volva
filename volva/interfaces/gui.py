"""
Volva XRD Analysis Tool - GUI Interface
Complete GUI application with SRM instrumental broadening correction
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from datetime import datetime
from typing import Dict, List, Optional, Any
import logging

from ..core import (
    VolvaConfig, Constants, AnalysisError, DataLoadError,
    PeakDetectionError, ValidationError, LMFIT_AVAILABLE
)
from ..data import XRDDataProcessor
from ..peaks import PeakDetector, PeakFitResult
from ..analysis import (
    ScherrerAnalyzer, MonshiScherrerAnalyzer,
    WilliamsonHallAnalyzer, HalderWagnerAnalyzer
)
from ..srm import (
    SRMProcessor, InstrumentalCorrector,
    get_srm, list_available_srms
)

logger = logging.getLogger(__name__)


class VolvaGUI:
    """Enhanced main GUI application with SRM instrumental correction support"""

    def __init__(self, config_file: str = "volva_config.json"):
        """Initialize the GUI application"""
        # Load configuration
        try:
            self.config = VolvaConfig.load(config_file)
        except Exception as e:
            logger.warning(f"Using default config due to error: {e}")
            self.config = VolvaConfig()

        self.config_file = config_file

        # Initialize main window
        self.root = tk.Tk()
        self.root.title("Volva - XRD Analysis Tool v4.0 Enhanced (SRM)")
        self.root.geometry("1600x1000")
        self.root.minsize(1200, 800)

        # Configure style
        self.setup_styles()

        # Data storage
        self.sample_data = {'two_theta': None, 'intensity': None, 'filename': None}
        self.srm_data = {'two_theta': None, 'intensity': None, 'filename': None}
        self.peaks = []
        self.srm_peaks = []
        self.analysis_results = {}

        # Current settings (needed before creating analyzers)
        self.current_wavelength_name = self.config.default_wavelength
        self.current_wavelength = Constants.WAVELENGTHS[self.current_wavelength_name]["ka1"]
        self.current_wavelength_info = Constants.WAVELENGTHS[self.current_wavelength_name]

        # Initialize analysis components FIRST (before SRM components that depend on them)
        self.peak_detector = PeakDetector(self.config)
        self.scherrer_analyzer = ScherrerAnalyzer(self.current_wavelength, self.config)
        self.monshi_analyzer = MonshiScherrerAnalyzer(self.current_wavelength_info, self.config)
        self.wh_analyzer = WilliamsonHallAnalyzer(self.current_wavelength, self.config)
        self.hw_analyzer = HalderWagnerAnalyzer(self.current_wavelength, self.config)

        # SRM components (now that peak_detector exists)
        self.srm_processor = SRMProcessor(self.peak_detector)
        self.instrumental_corrector = InstrumentalCorrector()
        self.srm_calibration = None
        self.corrected_peaks = []

        # GUI setup
        self.setup_gui()

        # Bind close event
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        logger.info("Volva GUI with SRM support initialized successfully")

    def setup_styles(self):
        """Setup GUI styles and themes"""
        style = ttk.Style()
        try:
            style.theme_use('clam')
        except:
            pass
        style.configure("Title.TLabel", font=('TkDefaultFont', 12, 'bold'))
        style.configure("Status.TLabel", font=('TkDefaultFont', 9))
        style.configure("Accent.TButton", foreground='blue')
        style.configure("SRM.TButton", foreground='darkgreen')

    def setup_gui(self):
        """Setup the complete GUI layout"""
        self.setup_menu()

        # Main container
        main_paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left panel for controls
        left_panel = ttk.Frame(main_paned, width=400)
        main_paned.add(left_panel, weight=0)

        # Right panel for plot and results
        right_paned = ttk.PanedWindow(main_paned, orient=tk.VERTICAL)
        main_paned.add(right_paned, weight=1)

        # Plot panel (top right)
        plot_panel = ttk.Frame(right_paned)
        right_paned.add(plot_panel, weight=2)

        # Results panel (bottom right)
        results_panel = ttk.Frame(right_paned)
        right_paned.add(results_panel, weight=1)

        # Setup individual panels
        self.setup_controls(left_panel)
        self.setup_plot(plot_panel)
        self.setup_results(results_panel)
        self.setup_status_bar()

    def setup_menu(self):
        """Setup menu bar"""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Load XRD Data...", command=self.load_data)
        file_menu.add_command(label="Load SRM Data...", command=self.load_srm_data)
        file_menu.add_separator()
        file_menu.add_command(label="Export Results...", command=self.export_results)
        file_menu.add_command(label="Exit", command=self.on_closing)

        # Analysis menu
        analysis_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)
        analysis_menu.add_command(label="Detect Peaks", command=self.detect_peaks)
        analysis_menu.add_command(label="Run All Methods", command=self.run_all_methods)

        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command=self.show_about)

    def setup_controls(self, parent):
        """Setup control panel"""
        notebook = ttk.Notebook(parent)
        notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # SRM tab (highest priority)
        srm_frame = ttk.Frame(notebook)
        notebook.add(srm_frame, text="SRM Data")
        self.setup_srm_tab(srm_frame)

        # Data tab
        data_frame = ttk.Frame(notebook)
        notebook.add(data_frame, text="Sample Data")
        self.setup_data_tab(data_frame)

        # Analysis tab
        analysis_frame = ttk.Frame(notebook)
        notebook.add(analysis_frame, text="Analysis")
        self.setup_analysis_tab(analysis_frame)

    def setup_srm_tab(self, parent):
        """Setup SRM (Standard Reference Material) tab"""
        # SRM file loading
        srm_file_frame = ttk.LabelFrame(parent, text="SRM Data Loading", padding=10)
        srm_file_frame.pack(fill=tk.X, pady=5)

        ttk.Button(srm_file_frame, text="Load SRM Data", command=self.load_srm_data,
                  style="SRM.TButton").pack(fill=tk.X, pady=2)

        self.srm_file_label = ttk.Label(srm_file_frame, text="No SRM data loaded",
                                       foreground="gray")
        self.srm_file_label.pack(fill=tk.X, pady=2)

        # SRM material selection
        material_frame = ttk.LabelFrame(parent, text="SRM Material", padding=10)
        material_frame.pack(fill=tk.X, pady=5)

        ttk.Label(material_frame, text="Reference Material:").pack(anchor=tk.W)
        self.srm_material_var = tk.StringVar()
        try:
            srm_materials = list_available_srms()
        except:
            srm_materials = ["LaB6", "Silicon", "Al2O3"]

        srm_combo = ttk.Combobox(material_frame, textvariable=self.srm_material_var,
                                values=srm_materials, state="readonly")
        srm_combo.pack(fill=tk.X, pady=2)

        ttk.Button(material_frame, text="Auto-Detect Material",
                  command=self.auto_detect_srm).pack(fill=tk.X, pady=2)

        # SRM processing
        process_frame = ttk.LabelFrame(parent, text="SRM Processing", padding=10)
        process_frame.pack(fill=tk.X, pady=5)

        ttk.Button(process_frame, text="Process SRM Peaks", command=self.process_srm,
                  style="SRM.TButton").pack(fill=tk.X, pady=2)

        self.srm_status_label = ttk.Label(process_frame, text="No SRM processing",
                                         foreground="gray")
        self.srm_status_label.pack(fill=tk.X, pady=2)

        # Instrumental calibration
        calib_frame = ttk.LabelFrame(parent, text="Instrumental Calibration", padding=10)
        calib_frame.pack(fill=tk.X, pady=5)

        ttk.Button(calib_frame, text="Create Calibration", command=self.create_srm_calibration,
                  style="SRM.TButton").pack(fill=tk.X, pady=2)

        self.calibration_status_label = ttk.Label(calib_frame, text="No calibration",
                                                 foreground="gray")
        self.calibration_status_label.pack(fill=tk.X, pady=2)

        # Peak correction
        correction_frame = ttk.LabelFrame(parent, text="Peak Corrections", padding=10)
        correction_frame.pack(fill=tk.BOTH, expand=True, pady=5)

        ttk.Button(correction_frame, text="Apply SRM Corrections",
                  command=self.apply_srm_corrections, style="Accent.TButton").pack(fill=tk.X, pady=5)

        # SRM info display
        self.srm_info_text = scrolledtext.ScrolledText(correction_frame, height=6,
                                                      wrap=tk.WORD, font=('Consolas', 9))
        self.srm_info_text.pack(fill=tk.BOTH, expand=True, pady=2)
        self.srm_info_text.insert(tk.END, "Load SRM data to begin instrumental calibration")
        self.srm_info_text.config(state=tk.DISABLED)

    def setup_data_tab(self, parent):
        """Setup data loading tab"""
        # File loading
        file_frame = ttk.LabelFrame(parent, text="Sample Data Loading", padding=10)
        file_frame.pack(fill=tk.X, pady=5)

        ttk.Button(file_frame, text="Load Sample Data", command=self.load_data,
                  style="Accent.TButton").pack(fill=tk.X, pady=2)

        self.file_label = ttk.Label(file_frame, text="No file loaded", foreground="gray")
        self.file_label.pack(fill=tk.X, pady=2)

        # Data information
        info_frame = ttk.LabelFrame(parent, text="Data Information", padding=10)
        info_frame.pack(fill=tk.X, pady=5)

        self.data_info_text = scrolledtext.ScrolledText(info_frame, height=4,
                                                       wrap=tk.WORD, font=('Consolas', 9))
        self.data_info_text.pack(fill=tk.BOTH, expand=True)
        self.data_info_text.insert(tk.END, "Load data to see information")
        self.data_info_text.config(state=tk.DISABLED)

        # Peak detection with enhanced controls
        peak_frame = ttk.LabelFrame(parent, text="Peak Detection", padding=10)
        peak_frame.pack(fill=tk.X, pady=5)

        # Background method
        ttk.Label(peak_frame, text="Background Method:").pack(anchor=tk.W)
        self.bg_method_var = tk.StringVar(value=self.config.default_bg_method)
        bg_combo = ttk.Combobox(peak_frame, textvariable=self.bg_method_var,
                               values=["polynomial", "rolling_minimum", "linear"],
                               state="readonly")
        bg_combo.pack(fill=tk.X, pady=2)

        # Prominence factor
        ttk.Label(peak_frame, text="Prominence Factor:").pack(anchor=tk.W, pady=(10, 0))
        self.prominence_var = tk.DoubleVar(value=self.config.default_prominence_factor)
        prominence_frame = ttk.Frame(peak_frame)
        prominence_frame.pack(fill=tk.X, pady=2)

        prominence_scale = ttk.Scale(prominence_frame, from_=0.01, to=0.5,
                                   variable=self.prominence_var, orient=tk.HORIZONTAL)
        prominence_scale.pack(side=tk.LEFT, fill=tk.X, expand=True)

        self.prominence_value = ttk.Label(prominence_frame, text="", width=8)
        self.prominence_value.pack(side=tk.RIGHT, padx=(5, 0))

        def update_prominence_label(*args):
            self.prominence_value.config(text=f"{self.prominence_var.get():.3f}")
        self.prominence_var.trace('w', update_prominence_label)
        update_prominence_label()

        # Minimum width
        ttk.Label(peak_frame, text="Minimum Width (°):").pack(anchor=tk.W, pady=(10, 0))
        self.min_width_var = tk.DoubleVar(value=self.config.default_min_width_deg)
        width_frame = ttk.Frame(peak_frame)
        width_frame.pack(fill=tk.X, pady=2)

        width_scale = ttk.Scale(width_frame, from_=0.01, to=0.5,
                              variable=self.min_width_var, orient=tk.HORIZONTAL)
        width_scale.pack(side=tk.LEFT, fill=tk.X, expand=True)

        self.width_value = ttk.Label(width_frame, text="", width=8)
        self.width_value.pack(side=tk.RIGHT, padx=(5, 0))

        def update_width_label(*args):
            self.width_value.config(text=f"{self.min_width_var.get():.3f}°")
        self.min_width_var.trace('w', update_width_label)
        update_width_label()

        # Advanced fitting checkbox
        self.advanced_fitting_var = tk.BooleanVar(value=self.config.use_lmfit)
        ttk.Checkbutton(peak_frame, text="Use Advanced Peak Fitting (lmfit)",
                       variable=self.advanced_fitting_var).pack(anchor=tk.W, pady=(10, 5))

        ttk.Button(peak_frame, text="Detect Sample Peaks", command=self.detect_peaks,
                  style="Accent.TButton").pack(fill=tk.X, pady=5)

        self.peaks_label = ttk.Label(peak_frame, text="No peaks detected", foreground="gray")
        self.peaks_label.pack(fill=tk.X, pady=2)

    def setup_analysis_tab(self, parent):
        """Setup analysis tab"""
        methods_frame = ttk.LabelFrame(parent, text="Analysis Methods", padding=10)
        methods_frame.pack(fill=tk.X, pady=5)

        ttk.Button(methods_frame, text="Scherrer Analysis",
                  command=lambda: self.run_analysis('scherrer')).pack(fill=tk.X, pady=2)
        ttk.Button(methods_frame, text="Monshi-Scherrer Analysis",
                  command=lambda: self.run_analysis('monshi')).pack(fill=tk.X, pady=2)
        ttk.Button(methods_frame, text="Williamson-Hall Analysis",
                  command=lambda: self.run_analysis('wh')).pack(fill=tk.X, pady=2)
        ttk.Button(methods_frame, text="Halder-Wagner Analysis",
                  command=lambda: self.run_analysis('hw')).pack(fill=tk.X, pady=2)

        ttk.Separator(methods_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        ttk.Button(methods_frame, text="Run All Methods",
                  command=self.run_all_methods, style="Accent.TButton").pack(fill=tk.X, pady=5)

    def setup_plot(self, parent):
        """Setup matplotlib plot"""
        plot_frame = ttk.Frame(parent)
        plot_frame.pack(fill=tk.BOTH, expand=True)

        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar_frame = ttk.Frame(plot_frame)
        toolbar_frame.pack(fill=tk.X)
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()

        self.plot_empty()

    def setup_results(self, parent):
        """Setup results display"""
        results_frame = ttk.Frame(parent)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.results_text = scrolledtext.ScrolledText(results_frame, height=10,
                                                     wrap=tk.WORD, font=('Consolas', 9))
        self.results_text.pack(fill=tk.BOTH, expand=True)

        button_frame = ttk.Frame(results_frame)
        button_frame.pack(fill=tk.X, pady=5)
        ttk.Button(button_frame, text="Clear Results", command=self.clear_results).pack(side=tk.LEFT)

    def setup_status_bar(self):
        """Setup status bar"""
        status_frame = ttk.Frame(self.root)
        status_frame.pack(side=tk.BOTTOM, fill=tk.X)

        self.status_var = tk.StringVar(value="Volva v4.0 Ready")
        status_bar = ttk.Label(status_frame, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=2, pady=2)

    def load_data(self):
        """Load XRD data"""
        file_path = filedialog.askopenfilename(
            title="Select XRD Data File",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )

        if not file_path:
            return

        try:
            two_theta, intensity, message = XRDDataProcessor.load_xrd_file(file_path)

            self.sample_data = {
                'two_theta': two_theta,
                'intensity': intensity,
                'filename': os.path.basename(file_path)
            }

            self.file_label.config(text=f"✅ {os.path.basename(file_path)}", foreground="darkgreen")
            self.plot_data()
            self.status_var.set(message)

            # Reset analysis state
            self.peaks = []
            self.analysis_results = {}
            self.clear_results()
            self.update_data_info()

        except Exception as e:
            error_msg = f"Data loading failed: {str(e)}"
            self.status_var.set(error_msg)
            messagebox.showerror("Error", error_msg)

    def load_srm_data(self):
        """Load SRM data"""
        file_path = filedialog.askopenfilename(
            title="Select SRM Data File",
            filetypes=[("Text files", "*.txt"), ("Data files", "*.dat *.xy"), ("All files", "*.*")]
        )

        if not file_path:
            return

        try:
            self.status_var.set("Loading SRM data...")
            self.root.update()

            # Use SRM processor to load data
            two_theta, intensity, metadata = self.srm_processor.load_srm_data(file_path)

            self.srm_data = {
                'two_theta': two_theta,
                'intensity': intensity,
                'filename': os.path.basename(file_path),
                'metadata': metadata
            }

            self.srm_file_label.config(text=f"✅ SRM: {os.path.basename(file_path)}",
                                      foreground="darkgreen")
            self.update_srm_info("SRM data loaded successfully")
            self.status_var.set("SRM data loaded successfully")

        except Exception as e:
            error_msg = f"SRM loading failed: {str(e)}"
            self.status_var.set(error_msg)
            messagebox.showerror("Error", error_msg)

    def auto_detect_srm(self):
        """Auto-detect SRM material"""
        if self.srm_data['two_theta'] is None:
            messagebox.showwarning("Warning", "Please load SRM data first")
            return

        try:
            self.status_var.set("Auto-detecting SRM material...")
            self.root.update()

            detected_material = self.srm_processor.identify_srm_material(
                self.srm_data['two_theta'],
                self.srm_data['intensity'],
                self.current_wavelength
            )

            self.srm_material_var.set(detected_material)
            self.update_srm_info(f"Auto-detected SRM material: {detected_material}")
            self.status_var.set(f"SRM material detected: {detected_material}")

        except Exception as e:
            error_msg = f"SRM auto-detection failed: {str(e)}"
            self.status_var.set(error_msg)
            messagebox.showerror("Error", error_msg)

    def process_srm(self):
        """Process SRM data"""
        if self.srm_data['two_theta'] is None:
            messagebox.showwarning("Warning", "Please load SRM data first")
            return

        if not self.srm_material_var.get():
            messagebox.showwarning("Warning", "Please select SRM material")
            return

        try:
            self.status_var.set("Processing SRM peaks...")
            self.root.update()

            self.srm_result = self.srm_processor.process_srm_data(
                self.srm_data['two_theta'],
                self.srm_data['intensity'],
                self.srm_material_var.get(),
                self.current_wavelength
            )

            matched_peaks = len(self.srm_result.matched_peaks)
            total_expected = matched_peaks + len(self.srm_result.missing_expected)

            self.srm_status_label.config(
                text=f"✅ {matched_peaks}/{total_expected} SRM peaks matched",
                foreground="darkgreen"
            )

            info_text = f"SRM Processing Results:\n"
            info_text += f"Material: {self.srm_result.srm_name}\n"
            info_text += f"Matched peaks: {matched_peaks}\n"
            info_text += f"Quality score: {self.srm_result.quality_score:.3f}\n"

            self.update_srm_info(info_text)
            self.status_var.set("SRM processing complete")

        except Exception as e:
            error_msg = f"SRM processing failed: {str(e)}"
            self.status_var.set(error_msg)
            messagebox.showerror("Error", error_msg)

    def create_srm_calibration(self):
        """Create SRM calibration"""
        if not hasattr(self, 'srm_result'):
            messagebox.showwarning("Warning", "Please process SRM data first")
            return

        try:
            self.status_var.set("Creating SRM calibration...")
            self.root.update()

            instrument_config = {
                'wavelength': self.current_wavelength,
                'radiation': self.current_wavelength_name
            }

            self.srm_calibration = self.instrumental_corrector.create_calibration(
                self.srm_result,
                instrument_config,
                function_name="Caglioti",
                auto_select_function=True
            )

            function_used = self.srm_calibration.fit_result.function_name
            r_squared = self.srm_calibration.fit_result.r_squared

            self.calibration_status_label.config(
                text=f"✅ {function_used} calibration (R² = {r_squared:.4f})",
                foreground="darkgreen"
            )

            calib_info = f"Instrumental Calibration Created:\n"
            calib_info += f"Function: {function_used}\n"
            calib_info += f"R²: {r_squared:.4f}\n"

            self.update_srm_info(calib_info)
            self.status_var.set("SRM calibration created successfully")

        except Exception as e:
            error_msg = f"SRM calibration failed: {str(e)}"
            self.status_var.set(error_msg)
            messagebox.showerror("Error", error_msg)

    def apply_srm_corrections(self):
        """Apply SRM corrections to sample peaks"""
        if not self.peaks:
            messagebox.showwarning("Warning", "Please detect sample peaks first")
            return

        if self.srm_calibration is None:
            messagebox.showwarning("Warning", "Please create SRM calibration first")
            return

        try:
            self.status_var.set("Applying SRM corrections...")
            self.root.update()

            # Convert peaks to PeakFitResult format
            peak_results = []
            for peak in self.peaks:
                peak_result = PeakFitResult(
                    fwhm=peak['FWHM_obs'],
                    r_squared=peak.get('r_squared', 0.85),
                    fit_info={'position': peak['2theta'], 'intensity': peak['intensity']},
                    center=peak['2theta'],
                    amplitude=peak['intensity']
                )
                peak_results.append(peak_result)

            # Apply corrections
            self.corrected_peaks = self.instrumental_corrector.apply_correction(
                peak_results,
                self.srm_calibration,
                method="deconvolution"
            )

            correction_info = f"SRM Corrections Applied:\n"
            correction_info += f"Peaks corrected: {len(self.corrected_peaks)}\n"

            self.update_srm_info(correction_info)
            self.status_var.set("SRM corrections applied successfully")

            # Update plot to show corrections
            if self.peaks:
                self.plot_data_with_peaks()

        except Exception as e:
            error_msg = f"SRM correction failed: {str(e)}"
            self.status_var.set(error_msg)
            messagebox.showerror("Error", error_msg)

    def update_srm_info(self, text: str):
        """Update SRM information display"""
        self.srm_info_text.config(state=tk.NORMAL)
        self.srm_info_text.delete(1.0, tk.END)
        self.srm_info_text.insert(tk.END, text)
        self.srm_info_text.config(state=tk.DISABLED)



    def detect_peaks(self):
        """Detect peaks with enhanced parameters"""
        if self.sample_data['two_theta'] is None:
            messagebox.showwarning("Warning", "Please load XRD data first")
            return

        try:
            self.status_var.set("Processing background subtraction...")
            self.root.update()

            # Background subtraction
            bg_method = self.bg_method_var.get()
            processed_intensity, background = XRDDataProcessor.subtract_background(
                self.sample_data['two_theta'],
                self.sample_data['intensity'],
                method=bg_method
            )

            self.status_var.set("Detecting peaks...")
            self.root.update()

            # Enhanced peak detection with parameters
            self.peaks = self.peak_detector.detect_peaks(
                self.sample_data['two_theta'],
                processed_intensity,
                prominence_factor=self.prominence_var.get(),
                min_width_deg=self.min_width_var.get(),
                use_advanced_fitting=self.advanced_fitting_var.get()
            )

            if self.peaks:
                # Quality filtering
                filtered_peaks = self.peak_detector.filter_peaks_by_quality(
                    self.peaks, self.config.min_r_squared)

                if len(filtered_peaks) < len(self.peaks):
                    quality_msg = f"Quality filtering: {len(filtered_peaks)}/{len(self.peaks)} peaks retained"
                    self.status_var.set(quality_msg)

                self.peaks = filtered_peaks
                self.peaks_label.config(text=f"✅ {len(self.peaks)} peaks detected",
                                        foreground="darkgreen")
                self.plot_data_with_peaks()
                self.status_var.set(f"Peak detection complete: {len(self.peaks)} peaks found")

                # Update data info with peak statistics
                self.update_data_info()
            else:
                self.peaks_label.config(text="❌ No peaks detected", foreground="red")
                messagebox.showinfo("Peak Detection",
                                    "No peaks detected. Try adjusting detection parameters.")
                self.status_var.set("Peak detection failed - no peaks found")

        except Exception as e:
            error_msg = f"Peak detection failed: {str(e)}"
            self.status_var.set(error_msg)
            messagebox.showerror("Error", error_msg)
            logger.error(f"Peak detection failed: {e}", exc_info=True)

    def update_data_info(self):
        """Update data information display"""
        if self.sample_data['two_theta'] is None:
            return

        two_theta = self.sample_data['two_theta']
        intensity = self.sample_data['intensity']

        info_text = f"Sample Data Information:\n"
        info_text += f"Filename: {self.sample_data['filename']}\n"
        info_text += f"Data points: {len(two_theta):,}\n"
        info_text += f"2θ range: {two_theta.min():.1f}° - {two_theta.max():.1f}°\n"
        info_text += f"Intensity range: {intensity.min():.0f} - {intensity.max():.0f} counts\n"
        info_text += f"Step size: {np.median(np.diff(two_theta)):.4f}°\n"

        if self.peaks:
            info_text += f"\nPeak Information:\n"
            info_text += f"Detected peaks: {len(self.peaks)}\n"
            avg_fwhm = np.mean([p['FWHM_obs'] for p in self.peaks])
            avg_quality = np.mean([p.get('r_squared', 0.5) for p in self.peaks])
            info_text += f"Average FWHM: {avg_fwhm:.3f}°\n"
            info_text += f"Average quality: {avg_quality:.3f}"

        self.data_info_text.config(state=tk.NORMAL)
        self.data_info_text.delete(1.0, tk.END)
        self.data_info_text.insert(tk.END, info_text)
        self.data_info_text.config(state=tk.DISABLED)

    def run_analysis(self, method_name: str):
        """Run specific analysis method with SRM support"""
        if not self.peaks:
            messagebox.showwarning("Warning", "Please detect peaks first")
            return

        try:
            self.status_var.set(f"Running {method_name} analysis...")
            self.root.update()

            # Use corrected peaks if available
            peaks_to_analyze = self.corrected_peaks if self.corrected_peaks else self.peaks

            if method_name == 'scherrer':
                result = self.scherrer_analyzer.analyze(peaks_to_analyze)
                display_name = "Scherrer Analysis"
            elif method_name == 'monshi':
                result = self.monshi_analyzer.analyze(peaks_to_analyze)
                display_name = "Monshi-Scherrer Analysis"
            elif method_name == 'wh':
                if len(peaks_to_analyze) < 3:
                    raise AnalysisError("Need at least 3 peaks for Williamson-Hall analysis")
                result = self.wh_analyzer.analyze(peaks_to_analyze)
                display_name = "Williamson-Hall Analysis"
            elif method_name == 'hw':
                if len(peaks_to_analyze) < 3:
                    raise AnalysisError("Need at least 3 peaks for Halder-Wagner analysis")
                result = self.hw_analyzer.analyze(peaks_to_analyze)
                display_name = "Halder-Wagner Analysis"
            else:
                raise ValueError(f"Unknown analysis method: {method_name}")

            if result['success']:
                self.analysis_results[method_name] = result
                self.display_results(display_name, result)
                self.status_var.set(f"{display_name} completed successfully")
            else:
                raise AnalysisError(result.get('message', 'Analysis failed'))

        except Exception as e:
            error_msg = f"Analysis failed: {str(e)}"
            self.status_var.set(error_msg)
            messagebox.showerror("Error", error_msg)
            logger.error(f"Analysis {method_name} failed: {e}", exc_info=True)

    def run_all_methods(self):
        """Run all available analysis methods"""
        if not self.peaks:
            messagebox.showwarning("Warning", "Please detect peaks first")
            return

        methods = ['scherrer', 'monshi']
        if len(self.peaks) >= 3:
            methods.extend(['wh', 'hw'])

        for method in methods:
            try:
                self.run_analysis(method)
            except Exception as e:
                logger.error(f"Method {method} failed: {e}")

    def plot_empty(self):
        """Plot empty state"""
        self.ax.clear()
        self.ax.set_xlabel('2theta (degrees)', fontsize=12)
        self.ax.set_ylabel('Intensity (counts)', fontsize=12)
        self.ax.set_title('Volva XRD Analysis Tool v4.0 Enhanced (SRM)', fontsize=14, fontweight='bold')

        empty_text = ('Load XRD data to begin analysis\n\n'
                      'SRM-Enhanced Features:\n'
                      '1. Load SRM Data (LaB₆ or Silicon)\n'
                      '2. Create instrumental calibration\n'
                      '3. Load sample data\n'
                      '4. Apply SRM corrections\n'
                      '5. Run analysis with corrected peaks')

        self.ax.text(0.5, 0.5, empty_text,
                     transform=self.ax.transAxes, ha='center', va='center',
                     fontsize=10, alpha=0.7, bbox=dict(boxstyle="round,pad=0.5",
                                                       facecolor="lightblue", alpha=0.1))
        self.ax.grid(True, alpha=0.3)
        self.canvas.draw()

    def plot_data(self):
        """Plot XRD data"""
        if self.sample_data['two_theta'] is None:
            self.plot_empty()
            return

        self.ax.clear()
        self.ax.plot(self.sample_data['two_theta'], self.sample_data['intensity'],
                     'b-', linewidth=1.5, label='Sample Data', alpha=0.8)

        self.ax.set_xlabel('2theta (degrees)', fontsize=12)
        self.ax.set_ylabel('Intensity (counts)', fontsize=12)
        self.ax.set_title(f'XRD Pattern: {self.sample_data["filename"]} ({self.current_wavelength_name})',
                          fontsize=13, fontweight='bold')
        self.ax.grid(True, alpha=0.3)
        self.ax.legend()

        # Add data statistics
        two_theta = self.sample_data['two_theta']
        stats_text = f"Points: {len(two_theta):,}\n"
        stats_text += f"Range: {two_theta.min():.1f} - {two_theta.max():.1f}°"

        if self.srm_calibration:
            stats_text += f"\nSRM: {self.srm_calibration.srm_name}"
        else:
            stats_text += f"\nSRM: Not applied"

        self.ax.text(0.02, 0.98, stats_text, transform=self.ax.transAxes,
                     va='top', fontsize=9,
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

        self.canvas.draw()

    def plot_data_with_peaks(self):
        """Plot XRD data with enhanced peak visualization"""
        if self.sample_data['two_theta'] is None or not self.peaks:
            return

        self.ax.clear()
        self.ax.plot(self.sample_data['two_theta'], self.sample_data['intensity'],
                     'b-', linewidth=1.5, label='Sample Data', alpha=0.8)

        # Plot peaks with quality-based coloring
        peak_positions = [p['2theta'] for p in self.peaks]
        peak_intensities = [p['intensity'] for p in self.peaks]
        peak_qualities = [p.get('r_squared', 0.5) for p in self.peaks]

        # Color-code peaks by quality
        colors = ['red' if q > 0.8 else 'orange' if q > 0.6 else 'yellow' for q in peak_qualities]

        self.ax.scatter(peak_positions, peak_intensities, c=colors, s=120,
                        marker='^', zorder=5, alpha=0.9, edgecolors='black', linewidth=1,
                        label='Sample Peaks')

        # Peak labels with FWHM info
        for i, (pos, intensity, fwhm, quality) in enumerate(zip(
                peak_positions, peak_intensities,
                [p['FWHM_obs'] for p in self.peaks], peak_qualities)):

            if self.corrected_peaks and i < len(self.corrected_peaks):
                corrected_fwhm = self.corrected_peaks[i].fwhm
                label = f'P{i + 1}\nObs: {fwhm:.3f}°\nCorr: {corrected_fwhm:.3f}°'
            else:
                label = f'P{i + 1}\n{fwhm:.3f}°'

            self.ax.annotate(label, (pos, intensity),
                             xytext=(5, 15), textcoords='offset points',
                             fontsize=7, fontweight='bold', ha='center',
                             bbox=dict(boxstyle="round,pad=0.2", facecolor=colors[i], alpha=0.3))

        self.ax.set_xlabel('2theta (degrees)', fontsize=12)
        self.ax.set_ylabel('Intensity (counts)', fontsize=12)

        # Enhanced title with SRM status
        title = f'Peak Detection: {self.sample_data["filename"]} ({len(self.peaks)} peaks)'
        if self.corrected_peaks:
            title += ' - SRM Corrected'

        self.ax.set_title(title, fontsize=13, fontweight='bold')
        self.ax.grid(True, alpha=0.3)

        # Enhanced legend
        from matplotlib.lines import Line2D
        legend_elements = [Line2D([0], [0], color='b', linewidth=1.5, label='Sample Data')]
        legend_elements.extend([
            Line2D([0], [0], marker='^', color='w', markerfacecolor='red', markersize=8,
                   label='High Quality (R² > 0.8)', markeredgecolor='black'),
            Line2D([0], [0], marker='^', color='w', markerfacecolor='orange', markersize=8,
                   label='Good Quality (R² > 0.6)', markeredgecolor='black'),
            Line2D([0], [0], marker='^', color='w', markerfacecolor='yellow', markersize=8,
                   label='Fair Quality (R² ≤ 0.6)', markeredgecolor='black')
        ])

        self.ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

        # Peak statistics with SRM info
        avg_fwhm = np.mean([p['FWHM_obs'] for p in self.peaks])
        avg_quality = np.mean(peak_qualities)
        stats_text = f"Peaks: {len(self.peaks)}\nAvg FWHM: {avg_fwhm:.3f}°\nAvg Quality: {avg_quality:.3f}"

        if self.corrected_peaks:
            avg_corrected = np.mean([p.fwhm for p in self.corrected_peaks])
            correction_factor = avg_corrected / avg_fwhm
            stats_text += f"\nSRM Corrected: {correction_factor:.2f}"

        self.ax.text(0.02, 0.98, stats_text, transform=self.ax.transAxes,
                     va='top', fontsize=9,
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.8))

        self.canvas.draw()

    def display_results(self, method_name: str, result: Dict):
        """Display enhanced analysis results"""
        text = f"\n{'=' * 60}\n{method_name.upper()}\n{'=' * 60}\n"
        text += f"Analysis Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        text += f"Sample: {self.sample_data.get('filename', 'Unknown')}\n"
        text += f"Wavelength: {result.get('wavelength', self.current_wavelength):.5f} nm\n"

        # SRM correction status
        if self.corrected_peaks:
            text += f"SRM Corrected: YES\n"
        else:
            text += f"SRM Corrected: NO\n"

        text += "\n"

        # Method-specific results
        if method_name.startswith("Scherrer") or method_name.startswith("Monshi"):
            text += f"Average Crystallite Size: {result['average_size']:.1f} ± {result['std_size']:.1f} nm\n"
            text += f"Median Crystallite Size: {result['median_size']:.1f} nm\n"
            text += f"Peaks Successfully Analyzed: {result['num_peaks']}\n"

            if method_name.startswith("Monshi"):
                text += f"Ka1 Wavelength: {result.get('ka1_wavelength', 0):.5f} nm\n"
                text += f"Ka2 Wavelength: {result.get('ka2_wavelength', 0):.5f} nm\n"

        elif method_name.startswith("Williamson") or method_name.startswith("Halder"):
            text += f"Crystallite Size: {result['crystallite_size_nm']:.1f}"
            if 'crystallite_size_uncertainty_nm' in result:
                text += f" ± {result['crystallite_size_uncertainty_nm']:.1f}"
            text += " nm\n"

            text += f"Microstrain: {result['microstrain_percent']:.5f}"
            if 'microstrain_uncertainty_percent' in result:
                text += f" ± {result['microstrain_uncertainty_percent']:.5f}"
            text += " %\n"

            text += f"Linear Fit R²: {result['r_squared']:.4f}\n"
            text += f"Peaks Used: {result['num_peaks']}\n"

        self.results_text.insert(tk.END, text + "\n")
        self.results_text.see(tk.END)

    def clear_results(self):
        """Clear results display"""
        self.results_text.delete(1.0, tk.END)
        self.analysis_results = {}

    def export_results(self):
        """Export analysis results"""
        if not self.analysis_results:
            messagebox.showwarning("Warning", "No results to export")
            return

        filename = filedialog.asksaveasfilename(
            title="Export Results",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w') as f:
                    f.write(self.results_text.get(1.0, tk.END))
                messagebox.showinfo("Export Success", f"Results exported to {filename}")
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export: {str(e)}")

    def show_about(self):
        """Show about dialog"""
        about_text = """Volva XRD Analysis Tool v4.0 Enhanced (SRM)

Advanced X-ray Diffraction analysis software with Standard Reference 
Material (SRM) instrumental broadening correction capability.

SRM FEATURES:
• LaB₆ and Silicon standard reference material support
• Angle-dependent instrumental profile correction
• Multiple instrumental functions (Caglioti, Polynomial, Pearson-VII)

ANALYSIS METHODS:
• Scherrer Analysis - Basic crystallite size calculation
• Monshi-Scherrer - Enhanced with Kα₂ correction  
• Williamson-Hall - Size and strain with isotropy validation
• Halder-Wagner - Anisotropic strain analysis

© 2025 Volva Development Team
Licensed under MIT License

Built with: Python, NumPy, SciPy, Matplotlib, Tkinter, lmfit"""

        messagebox.showinfo("About Volva SRM", about_text)

    def on_closing(self):
        """Handle application closing"""
        try:
            self.config.save(self.config_file)
            logger.info("Configuration saved on exit")
        except Exception as e:
            logger.warning(f"Failed to save configuration on exit: {e}")
        self.root.destroy()

    def run(self):
        """Start the GUI application"""
        try:
            self.root.update_idletasks()
            width = self.root.winfo_width()
            height = self.root.winfo_height()
            x = (self.root.winfo_screenwidth() // 2) - (width // 2)
            y = (self.root.winfo_screenheight() // 2) - (height // 2)
            self.root.geometry(f"{width}x{height}+{x}+{y}")

            logger.info("Starting Volva SRM GUI main loop")
            self.root.mainloop()

        except Exception as e:
            logger.error(f"GUI application error: {e}", exc_info=True)
            messagebox.showerror("Application Error", f"A critical error occurred: {str(e)}")

# Module-level launch function - NO INDENTATION (at module level)
def launch_gui(config_file: str = "volva_config.json"):
    """Launch the Volva SRM-enhanced GUI application"""
    try:
        app = VolvaGUI(config_file)
        app.run()
    except Exception as e:
        print(f"GUI launch failed: {e}")
        logger.error(f"GUI launch failed: {e}", exc_info=True)
        print("Make sure tkinter is available in your Python installation")
        if LMFIT_AVAILABLE:
            print("lmfit is available - advanced peak fitting will work")
        else:
            print("Consider installing lmfit for advanced peak fitting: pip install lmfit")