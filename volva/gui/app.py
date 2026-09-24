"""
Völva desktop application (PySide6).

Run from Anaconda Prompt:   python -m volva
"""

from __future__ import annotations

import html
import os
import sys
import traceback
import warnings

import numpy as np
from PySide6.QtCore import QObject, QSettings, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QAction, QFontDatabase, QIcon, QKeySequence
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow,
                               QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QScrollArea, QSpinBox,
                               QSplitter, QTableWidget, QTableWidgetItem, QTabWidget, QTextBrowser, QVBoxLayout,
                               QWidget)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from .. import __version__  # noqa: E402
from .. import constants as C  # noqa: E402
from .. import plots  # noqa: E402
from ..instrument import InstrumentProfile, calibrate_from_standard  # noqa: E402
from ..io import LoadError, load_pattern  # noqa: E402
from ..pipeline import AnalysisSettings, rerun_methods, resolve_wavelength, run_analysis  # noqa: E402
from ..methods import _bound  # noqa: E402
from ..report import text_report, write_bundle  # noqa: E402

LEVEL_STYLE = {
    "critical": ("CRITICAL", "#b00020"),
    "warning": ("WARNING", "#8a5a00"),
    "info": ("NOTE", "#1f4e79"),
    "ok": ("OK", "#1b5e20"),
}


# ---------------------------------------------------------------------------
class Worker(QObject):
    finished = Signal(object)
    failed = Signal(str)
    progress = Signal(str)

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn, self.args, self.kwargs = fn, args, kwargs

    def run(self):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                out = self.fn(*self.args, **self.kwargs)
            self.finished.emit(out)
        except Exception:
            self.failed.emit(traceback.format_exc())


class PlotTab(QWidget):
    def __init__(self, draw, parent=None, extra=None):
        super().__init__(parent)
        self.draw_fn = draw
        self.fig = Figure(figsize=(9, 5.5), facecolor="white")
        self.canvas = FigureCanvasQTAgg(self.fig)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        top = QHBoxLayout()
        top.addWidget(NavigationToolbar2QT(self.canvas, self))
        if extra is not None:
            top.addWidget(extra)
        top.addStretch()
        lay.addLayout(top)
        lay.addWidget(self.canvas, 1)

    def redraw(self):
        try:
            self.draw_fn(self.fig)
        except Exception as e:  # never let a plot crash the app
            plots._empty(self.fig, f"Plot error: {e}")
        self.canvas.draw_idle()


# ---------------------------------------------------------------------------
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"Völva {__version__} — crystallite size and strain")
        self.resize(1500, 950)
        self.qs = QSettings("UniversityOfSheffield", "Volva")
        self.pattern = None
        self.profile = None
        self.profile_path = ""
        self.result = None
        self.peak_index = 0
        self._thread = None
        self._table_updating = False
        self._rerun_timer = QTimer(self, singleShot=True, interval=250, timeout=self._rerun)

        self._build_menu()
        self._build_ui()
        self._autoload_profile()
        self._refresh_all()

    # ------------------------------------------------------------------ UI
    def _build_menu(self):
        m = self.menuBar().addMenu("&File")
        for text, slot, key in [
            ("&Open sample pattern…", self.open_pattern, QKeySequence.Open),
            ("&Calibrate from standard…", self.calibrate, None),
            ("&Load instrument profile…", self.load_profile, None),
            ("&Save instrument profile…", self.save_profile, None),
            (None, None, None),
            ("&Export results…", self.export_results, QKeySequence("Ctrl+E")),
            ("Save current &figure…", self.save_figure, QKeySequence.Save),
            (None, None, None),
            ("&Quit", self.close, QKeySequence.Quit),
        ]:
            if text is None:
                m.addSeparator()
                continue
            a = QAction(text, self)
            a.triggered.connect(slot)
            if key is not None:
                a.setShortcut(key)
            m.addAction(a)
        h = self.menuBar().addMenu("&Help")
        a = QAction("About Völva", self)
        a.triggered.connect(self.about)
        h.addAction(a)

    def _group(self, title):
        g = QGroupBox(title)
        f = QFormLayout(g)
        f.setLabelAlignment(Qt.AlignRight)
        return g, f

    def _build_ui(self):
        split = QSplitter(Qt.Horizontal)
        self.setCentralWidget(split)

        # ---------------- left panel
        left = QWidget()
        lv = QVBoxLayout(left)

        g, f = self._group("1. Instrument profile")
        self.lbl_profile = QLabel("None loaded")
        self.lbl_profile.setWordWrap(True)
        f.addRow(self.lbl_profile)
        self.cmb_standard = QComboBox()
        self.cmb_standard.addItems(list(C.STANDARDS))
        f.addRow("Standard:", self.cmb_standard)
        self.ed_note = QLineEdit()
        self.ed_note.setPlaceholderText("e.g. Aeris, Cu, 1/4° divergence slit")
        f.addRow("Optics note:", self.ed_note)
        row = QHBoxLayout()
        b1 = QPushButton("Calibrate…")
        b1.setToolTip("Fit a line-profile standard (e.g. NIST SRM 660c LaB₆) measured with the same optics")
        b1.clicked.connect(self.calibrate)
        b2 = QPushButton("Load…")
        b2.clicked.connect(self.load_profile)
        b3 = QPushButton("Save…")
        b3.clicked.connect(self.save_profile)
        for b in (b1, b2, b3):
            row.addWidget(b)
        f.addRow(row)
        lv.addWidget(g)

        g, f = self._group("2. Sample")
        self.lbl_sample = QLabel("No pattern loaded")
        self.lbl_sample.setWordWrap(True)
        f.addRow(self.lbl_sample)
        b = QPushButton("Open pattern…")
        b.clicked.connect(self.open_pattern)
        f.addRow(b)
        rr = QHBoxLayout()
        self.chk_range = QCheckBox("Limit 2θ to")
        self.sp_tmin = QDoubleSpinBox(decimals=2, minimum=1, maximum=179, value=15)
        self.sp_tmax = QDoubleSpinBox(decimals=2, minimum=1, maximum=179, value=140)
        rr.addWidget(self.chk_range)
        rr.addWidget(self.sp_tmin)
        rr.addWidget(QLabel("–"))
        rr.addWidget(self.sp_tmax)
        f.addRow(rr)
        lv.addWidget(g)

        g, f = self._group("3. Settings")
        self.cmb_wl = QComboBox()
        self.cmb_wl.addItems(list(C.WAVELENGTHS) + ["From file", "Custom"])
        self.cmb_wl.currentTextChanged.connect(self._wl_changed)
        f.addRow("Radiation:", self.cmb_wl)
        self.sp_ka1 = QDoubleSpinBox(decimals=7, minimum=0.01, maximum=1.0, singleStep=1e-4, value=0.1540562)
        self.sp_ka2 = QDoubleSpinBox(decimals=7, minimum=0.01, maximum=1.0, singleStep=1e-4, value=0.1544390)
        self.sp_ka1.setEnabled(False)
        self.sp_ka2.setEnabled(False)
        f.addRow("Kα1 / Kα2 (nm):", self._hbox(self.sp_ka1, self.sp_ka2))
        self.chk_doublet = QCheckBox("Fit Kα2 component")
        self.chk_doublet.setChecked(True)
        f.addRow(self.chk_doublet)
        self.cmb_size = QComboBox()
        self.cmb_size.addItems(list(C.SIZE_CONVENTIONS) + ["Custom K"])
        self.cmb_size.currentTextChanged.connect(self._k_changed)
        f.addRow("Size reported as:", self.cmb_size)
        self.sp_k = QDoubleSpinBox(decimals=3, minimum=0.5, maximum=3.0, singleStep=0.01, value=1.0)
        self.sp_k.setEnabled(False)
        self.sp_k.valueChanged.connect(self._k_changed)
        f.addRow("Custom K:", self.sp_k)
        self.chk_bg_auto = QCheckBox("Automatic")
        self.chk_bg_auto.setChecked(True)
        self.sp_bg = QDoubleSpinBox(decimals=2, minimum=0.2, maximum=15, singleStep=0.1, value=2.0, suffix="°")
        self.sp_bg.setEnabled(False)
        self.chk_bg_auto.toggled.connect(lambda on: self.sp_bg.setEnabled(not on))
        f.addRow("Background half-width:", self._hbox(self.chk_bg_auto, self.sp_bg))
        self.sp_prom = QDoubleSpinBox(decimals=1, minimum=2, maximum=100, value=C.DEFAULT_PROMINENCE_SIGMA, suffix=" σ")
        self.sp_prom.setToolTip("A feature is a peak candidate if its net area above a local baseline is at least "
                                "this many standard deviations of the counting noise (5 is a good default)")
        f.addRow("Detection threshold:", self.sp_prom)
        self.sp_r2 = QDoubleSpinBox(decimals=3, minimum=0.5, maximum=0.999, singleStep=0.01, value=C.DEFAULT_MIN_R2)
        f.addRow("Minimum fit R²:", self.sp_r2)
        self.sp_ratio = QDoubleSpinBox(decimals=2, minimum=1.05, maximum=3, singleStep=0.05,
                                       value=C.DEFAULT_RESOLUTION_RATIO)
        self.sp_ratio.setToolTip("β_obs/β_inst below this: reflection is near the resolution limit. If most are, only a lower bound is reported")
        f.addRow("Resolution limit:", self.sp_ratio)
        lv.addWidget(g)

        self.btn_run = QPushButton("Analyse")
        self.btn_run.setMinimumHeight(38)
        self.btn_run.setStyleSheet("QPushButton{font-weight:bold;background:#0072B2;color:white;border-radius:4px}"
                                   "QPushButton:disabled{background:#9ab}")
        self.btn_run.clicked.connect(self.analyse)
        lv.addWidget(self.btn_run)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        lv.addWidget(self.progress)

        g = QGroupBox("Peaks (untick to exclude; type hkl labels)")
        gl = QVBoxLayout(g)
        self.table = QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels(["Use", "2θ (°)", "hkl", "FWHM (°)", "η", "β/βᵢ", "Scherrer (nm)", "Flags"])
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemChanged.connect(self._table_changed)
        self.table.cellDoubleClicked.connect(self._table_dbl)
        gl.addWidget(self.table)
        row = QHBoxLayout()
        for text, fn in [("All", lambda: self._set_all(True)), ("None", lambda: self._set_all(False)),
                         ("Only clean", self._only_clean)]:
            b = QPushButton(text)
            b.clicked.connect(fn)
            row.addWidget(b)
        gl.addLayout(row)
        lv.addWidget(g, 1)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(left)
        scroll.setMinimumWidth(480)
        split.addWidget(scroll)

        # ---------------- right: tabs
        self.tabs = QTabWidget()
        split.addWidget(self.tabs)
        split.setStretchFactor(1, 1)
        split.setSizes([500, 1000])

        self.txt_result = QTextBrowser()
        self.tabs.addTab(self.txt_result, "Result")
        self.plot_tabs = {}

        def add(name, fn, extra=None):
            t = PlotTab(fn, extra=extra)
            self.plot_tabs[name] = t
            self.tabs.addTab(t, name)

        add("Pattern", lambda fg: plots.plot_pattern(fg, self.result) if self.result else
            (plots.plot_raw(fg, self.pattern) if self.pattern else plots._empty(fg, "Open a pattern")))
        nav = QWidget()
        nl = QHBoxLayout(nav)
        nl.setContentsMargins(0, 0, 0, 0)
        bp = QPushButton("◀")
        bp.clicked.connect(lambda: self._step_peak(-1))
        bn = QPushButton("▶")
        bn.clicked.connect(lambda: self._step_peak(1))
        self.lbl_peak = QLabel("")
        for w in (bp, self.lbl_peak, bn):
            nl.addWidget(w)
        add("Peak", lambda fg: plots.plot_peak(fg, self.result, self.peak_index), extra=nav)
        add("Gallery", lambda fg: plots.plot_gallery(fg, self.result))
        add("Instrument", lambda fg: plots.plot_instrument(fg, self.profile, self.result))
        add("Size by reflection", lambda fg: plots.plot_sizes(fg, self.result))
        add("Williamson-Hall", lambda fg: plots.plot_linear(fg, self.result, "wh"))
        add("Voigt W-H", lambda fg: plots.plot_linear(fg, self.result, "voigt_wh"))
        add("Halder-Wagner", lambda fg: plots.plot_linear(fg, self.result, "hw"))
        add("Monshi", lambda fg: plots.plot_linear(fg, self.result, "monshi"))
        add("Comparison", lambda fg: plots.plot_comparison(fg, self.result))
        mono = QFontDatabase.systemFont(QFontDatabase.FixedFont)
        self.txt_report = QPlainTextEdit(readOnly=True)
        self.txt_report.setFont(mono)
        self.tabs.addTab(self.txt_report, "Report")
        self.txt_log = QPlainTextEdit(readOnly=True)
        self.txt_log.setFont(mono)
        self.tabs.addTab(self.txt_log, "Log")
        self.tabs.currentChanged.connect(self._tab_changed)
        self.statusBar().showMessage("Ready")

    @staticmethod
    def _hbox(*ws):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        for x in ws:
            h.addWidget(x)
        return w

    # ------------------------------------------------------------------ helpers
    def log(self, text):
        self.txt_log.appendPlainText(text)

    def settings(self) -> AnalysisSettings:
        s = AnalysisSettings(
            wavelength=self.cmb_wl.currentText(),
            custom_ka1=self.sp_ka1.value(), custom_ka2=self.sp_ka2.value(), custom_ratio=0.5,
            doublet=self.chk_doublet.isChecked(),
            bg_half_width=None if self.chk_bg_auto.isChecked() else self.sp_bg.value(),
            prominence_sigma=self.sp_prom.value(), min_r2=self.sp_r2.value(),
            resolution_ratio=self.sp_ratio.value(),
            two_theta_min=self.sp_tmin.value() if self.chk_range.isChecked() else None,
            two_theta_max=self.sp_tmax.value() if self.chk_range.isChecked() else None)
        if self.cmb_size.currentText() == "Custom K":
            s.k_custom = self.sp_k.value()
        else:
            s.size_convention = self.cmb_size.currentText()
        return s

    def _busy(self, on, msg=""):
        self.btn_run.setEnabled(not on)
        self.progress.setVisible(on)
        if msg:
            self.statusBar().showMessage(msg)

    def _start(self, fn, done, *args, **kw):
        if self._thread is not None:
            return
        self._thread = QThread(self)
        self._worker = Worker(fn, *args, **kw)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(done)
        self._worker.failed.connect(self._failed)
        self._worker.finished.connect(self._cleanup)
        self._worker.failed.connect(self._cleanup)
        self._thread.start()

    def _cleanup(self, *_):
        if self._thread:
            self._thread.quit()
            self._thread.wait()
        self._thread = None
        self._busy(False)

    def _failed(self, tb):
        self.log(tb)
        QMessageBox.critical(self, "Völva", "The operation failed. Details are in the Log tab.\n\n"
                             + tb.strip().splitlines()[-1])

    def _dir(self):
        return self.qs.value("last_dir", os.path.expanduser("~"))

    # ------------------------------------------------------------------ actions
    def open_pattern(self):
        fp, _ = QFileDialog.getOpenFileName(self, "Open sample pattern", self._dir(),
                                            "Diffraction data (*.xy *.xye *.txt *.dat *.csv *.asc *.xrdml);;All files (*)")
        if not fp:
            return
        self.qs.setValue("last_dir", os.path.dirname(fp))
        try:
            self.pattern = load_pattern(fp)
        except LoadError as e:
            QMessageBox.warning(self, "Could not load", str(e))
            return
        p = self.pattern
        self.result = None
        txt = f"{p.name}\n{len(p.two_theta)} points, {p.two_theta.min():.2f}–{p.two_theta.max():.2f}°, step {p.step:.4f}°"
        if p.wavelength:
            self.cmb_wl.setCurrentText("From file")
            txt += f"\nWavelength from file: Kα1 {p.wavelength['ka1']:.6f} nm"
        self.lbl_sample.setText(txt)
        for k, v in p.meta.items():
            if k.startswith("warning"):
                self.log("⚠ " + v)
        self.log(f"Loaded {p.path}")
        self._refresh_all()
        self.tabs.setCurrentWidget(self.plot_tabs["Pattern"])

    def calibrate(self):
        fp, _ = QFileDialog.getOpenFileName(self, "Open line-profile standard pattern", self._dir(),
                                            "Diffraction data (*.xy *.xye *.txt *.dat *.csv *.asc *.xrdml);;All files (*)")
        if not fp:
            return
        try:
            pat = load_pattern(fp)
        except LoadError as e:
            QMessageBox.warning(self, "Could not load", str(e))
            return
        s = self.settings()
        wl = resolve_wavelength(s, pat)
        self._busy(True, "Calibrating instrument profile…")
        self._start(calibrate_from_standard, self._calibrated, pat, wl, doublet=s.doublet,
                    standard=self.cmb_standard.currentText(), note=self.ed_note.text())

    def _calibrated(self, out):
        prof, log, fits, bg = out
        self.log("\n".join(log))
        if prof is None:
            self.statusBar().showMessage("Calibration failed — see the Log tab.")
            QMessageBox.warning(self, "Calibration failed", "Too few usable standard peaks. See the Log tab.")
            return
        self.profile = prof
        n_used = sum(1 for q in prof.peaks if q.get("used", True))
        self.statusBar().showMessage(
            f"Instrument profile calibrated from {prof.source_file} ({n_used} peaks"
            f"{', with warnings' if prof.warnings else ''}). Not yet saved.")
        self.profile_path = ""
        if self.result is not None:
            self.result = None
        self._show_profile()
        self._refresh_all()
        self.tabs.setCurrentWidget(self.plot_tabs["Instrument"])
        msg = prof.summary()
        if prof.warnings:
            QMessageBox.warning(self, "Calibration warnings", msg)
        if QMessageBox.question(self, "Save profile", "Save this instrument profile for future sessions?") \
                == QMessageBox.Yes:
            self.save_profile()

    def load_profile(self):
        fp, _ = QFileDialog.getOpenFileName(self, "Load instrument profile", self._dir(), "Völva profile (*.json)")
        if fp:
            self._load_profile_file(fp)

    def _load_profile_file(self, fp, quiet=False):
        try:
            self.profile = InstrumentProfile.load(fp)
        except Exception as e:
            if not quiet:
                QMessageBox.warning(self, "Could not load profile", str(e))
            return
        self.profile_path = fp
        self.qs.setValue("last_profile", fp)
        self.log(f"Instrument profile loaded: {fp}\n{self.profile.summary()}")
        self.result = None
        self._show_profile()
        self._refresh_all()
        self.statusBar().showMessage(self._profile_status(fp, "Loaded instrument profile"))

    def _profile_status(self, fp, prefix):
        msg = f"{prefix}: {os.path.basename(fp)}"
        if self.profile is not None and self.profile.profile_model != "total":
            msg += "  —  made with the Völva 4.1.0 model; please recalibrate from the standard"
        return msg

    def _autoload_profile(self):
        fp = self.qs.value("last_profile", "")
        if fp and os.path.isfile(fp):
            self._load_profile_file(fp, quiet=True)
            if self.profile is not None:
                self.statusBar().showMessage(self._profile_status(fp, "Loaded last instrument profile"))

    def save_profile(self):
        if self.profile is None:
            QMessageBox.information(self, "Völva", "No instrument profile to save.")
            return
        base = os.path.splitext(self.profile.source_file)[0] + "_profile.json"
        fp, _ = QFileDialog.getSaveFileName(self, "Save instrument profile", os.path.join(self._dir(), base),
                                            "Völva profile (*.json)")
        if not fp:
            return
        if self.ed_note.text() and not self.profile.note:
            self.profile.note = self.ed_note.text()
        self.profile.save(fp)
        self.profile_path = fp
        self.qs.setValue("last_profile", fp)
        self._show_profile()
        self.log(f"Saved instrument profile: {fp}")
        self.statusBar().showMessage(f"Instrument profile saved: {os.path.basename(fp)}")

    def _show_profile(self):
        p = self.profile
        if p is None:
            self.lbl_profile.setText("None loaded — sizes will be lower bounds")
            return
        t = f"{p.name}\n{p.standard}, {p.two_theta_range[0]:.0f}–{p.two_theta_range[1]:.0f}°, " \
            f"{len(p.peaks)} peaks, created {p.created[:10]}"
        if p.note:
            t += f"\n{p.note}"
        if p.warnings:
            t += f"\n⚠ {len(p.warnings)} warning(s) — see Log"
        self.lbl_profile.setText(t)

    def analyse(self):
        if self.pattern is None:
            QMessageBox.information(self, "Völva", "Open a sample pattern first.")
            return
        if self.profile is None:
            r = QMessageBox.question(
                self, "No instrument profile",
                "No instrument profile is loaded, so instrumental broadening cannot be removed and every size "
                "will be a lower bound.\n\nAnalyse anyway?")
            if r != QMessageBox.Yes:
                return
        s = self.settings()
        self._busy(True, "Analysing…")
        self._start(run_analysis, self._analysed, self.pattern, self.profile, s)

    def _analysed(self, res):
        self.result = res
        self.peak_index = 0
        self.log("\n".join(res.log))
        for w in res.warnings:
            self.log("⚠ " + w)
        self._fill_table()
        self._refresh_all()
        self.tabs.setCurrentWidget(self.txt_result)
        self.statusBar().showMessage(res.interpretation.headline if res.interpretation else "Done")

    def _rerun(self):
        if self.result is None:
            return
        use, hkl = [], []
        for i in range(self.table.rowCount()):
            use.append(self.table.item(i, 0).checkState() == Qt.Checked)
            hkl.append(self.table.item(i, 2).text().strip())
        s = self.result.settings
        s2 = self.settings()
        s.size_convention, s.k_custom = s2.size_convention, s2.k_custom
        rerun_methods(self.result, use=use, hkl=hkl, settings=s)
        self._fill_table(keep_state=True)
        self._refresh_all()
        if self.result.interpretation:
            self.statusBar().showMessage(self.result.interpretation.headline)

    def export_results(self):
        if self.result is None:
            QMessageBox.information(self, "Völva", "Nothing to export yet.")
            return
        d = QFileDialog.getExistingDirectory(self, "Choose a folder for the results", self._dir())
        if not d:
            return
        paths = write_bundle(self.result, d)
        self.log("Exported:\n  " + "\n  ".join(paths))
        QMessageBox.information(self, "Exported", f"{len(paths)} files written to\n{d}")

    def save_figure(self):
        w = self.tabs.currentWidget()
        if not isinstance(w, PlotTab):
            QMessageBox.information(self, "Völva", "Select a plot tab first.")
            return
        base = os.path.splitext(self.pattern.name if self.pattern else "volva")[0]
        name = self.tabs.tabText(self.tabs.currentIndex()).replace(" ", "_")
        fp, _ = QFileDialog.getSaveFileName(self, "Save figure", os.path.join(self._dir(), f"{base}_{name}.png"),
                                            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg);;TIFF (*.tiff)")
        if fp:
            w.fig.savefig(fp, dpi=600, facecolor="white",
                          metadata={"Software": f"Völva {__version__}"} if fp.lower().endswith((".png", ".pdf", ".svg")) else None)
            self.log(f"Saved figure: {fp}")

    def about(self):
        QMessageBox.about(self, "About Völva",
                          f"<b>Völva {__version__}</b><br>XRD crystallite size and microstrain analysis.<br><br>"
                          "Voigt-based instrumental correction, integral breadths, six complementary methods and "
                          "an explicit comparison of their assumptions.<br><br>University of Sheffield.")

    # ------------------------------------------------------------------ table
    def _fill_table(self, keep_state=False):
        res = self.result
        self._table_updating = True
        self.table.setRowCount(len(res.peaks) if res else 0)
        if res:
            K, lam = res.settings.K, res.wavelength["ka1"]
            for i, p in enumerate(res.peaks):
                it = QTableWidgetItem()
                it.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
                it.setCheckState(Qt.Checked if p.use else Qt.Unchecked)
                self.table.setItem(i, 0, it)
                vals = [f"{p.two_theta:.3f}", p.hkl, f"{p.fit.fwhm:.4f}", f"{p.fit.eta:.2f}", f"{p.ratio:.2f}"]
                if p.resolution_limited:
                    b = _bound(p, K, lam)
                    vals.append(f"> {b:.0f}" if np.isfinite(b) else "RL")
                else:
                    vals.append(f"{K * lam / (p.beta_s * np.cos(p.theta)):.1f}" if p.beta_s > 0 else "—")
                vals.append(", ".join(x.replace("_", " ") for x in p.flags + p.fit.flags))
                for j, v in enumerate(vals, start=1):
                    item = QTableWidgetItem(v)
                    if j == 2:
                        item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsEditable | Qt.ItemIsSelectable)
                    else:
                        item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                    if j in (1, 3, 4, 5, 6):
                        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                    self.table.setItem(i, j, item)
        self._table_updating = False

    def _table_changed(self, item):
        if self._table_updating:
            return
        if item.column() in (0, 2):
            self._rerun_timer.start()

    def _table_dbl(self, row, col):
        if col != 2:
            self.peak_index = row
            self.tabs.setCurrentWidget(self.plot_tabs["Peak"])
            self._refresh_all()

    def _set_all(self, on):
        self._table_updating = True
        for i in range(self.table.rowCount()):
            self.table.item(i, 0).setCheckState(Qt.Checked if on else Qt.Unchecked)
        self._table_updating = False
        self._rerun()

    def _only_clean(self):
        if not self.result:
            return
        self._table_updating = True
        for i, p in enumerate(self.result.peaks):
            bad = p.resolution_limited or any(f in p.fit.flags for f in ("overlapped", "eta_at_bound", "no_uncertainty"))
            self.table.item(i, 0).setCheckState(Qt.Unchecked if bad else Qt.Checked)
        self._table_updating = False
        self._rerun()

    # ------------------------------------------------------------------ refresh
    def _wl_changed(self, text):
        custom = text == "Custom"
        self.sp_ka1.setEnabled(custom)
        self.sp_ka2.setEnabled(custom)
        if text in C.WAVELENGTHS:
            self.sp_ka1.setValue(C.WAVELENGTHS[text]["ka1"])
            self.sp_ka2.setValue(C.WAVELENGTHS[text]["ka2"])
        if self.result is not None:
            self.statusBar().showMessage("Radiation changed — press Analyse to update.")

    def _k_changed(self, *_):
        self.sp_k.setEnabled(self.cmb_size.currentText() == "Custom K")
        if self.result is not None:
            self._rerun_timer.start()

    def _step_peak(self, d):
        if self.result and self.result.peaks:
            self.peak_index = int(np.clip(self.peak_index + d, 0, len(self.result.peaks) - 1))
            self._redraw_current()

    def _tab_changed(self, _):
        self._redraw_current()

    def _redraw_current(self):
        w = self.tabs.currentWidget()
        if isinstance(w, PlotTab):
            w.redraw()
        if self.result and self.result.peaks:
            self.lbl_peak.setText(f"{self.peak_index + 1}/{len(self.result.peaks)}")

    def _refresh_all(self):
        self._render_result()
        if self.result:
            self.txt_report.setPlainText(text_report(self.result))
        else:
            self.txt_report.setPlainText("")
        self._redraw_current()

    def _render_result(self):
        res = self.result
        if res is None:
            self.txt_result.setHtml(
                "<h2>Völva " + __version__ + "</h2>"
                "<ol><li><b>Instrument profile</b>: calibrate from a LaB₆ (SRM 660c) pattern measured with the "
                "same optics, or load a saved profile. " + ("<i>Loaded.</i>" if self.profile else "") + "</li>"
                "<li><b>Open</b> the sample pattern (.xy, .xye, .csv, .txt or .xrdml).</li>"
                "<li>Press <b>Analyse</b>. Check the Pattern and Peak tabs, untick doubtful peaks, add hkl labels.</li>"
                "<li>Read the <b>Result</b> tab: it compares the methods and says which assumptions hold.</li>"
                "<li><b>File ▸ Export results</b> writes the report, CSV tables, figures and a session file "
                "recording every setting and file checksum.</li></ol>")
            return
        it = res.interpretation
        H = [f"<h2 style='margin-bottom:4px'>{html.escape(res.pattern.name)}</h2>",
             f"<div style='border:1px solid #0072B2;border-radius:6px;padding:10px;background:#f2f8fc;"
             f"font-size:14px'><b>{html.escape(it.headline)}</b></div>"]
        if res.warnings:
            H.append("<h3>Setup warnings</h3><ul>")
            H += [f"<li style='color:#8a5a00'><b>WARNING</b> {html.escape(w)}</li>" for w in res.warnings]
            H.append("</ul>")
        H.append("<h3>What the methods say together</h3><table cellspacing='0' cellpadding='4'>")
        for f in it.findings:
            tag, col = LEVEL_STYLE.get(f.level, ("NOTE", "#333"))
            H.append(f"<tr><td valign='top' style='color:{col};font-weight:bold'>{tag}</td>"
                     f"<td valign='top'><b>{html.escape(f.topic)}</b></td><td>{html.escape(f.text)}</td></tr>")
        H.append("</table><h3>Methods</h3>"
                 "<table border='1' cellspacing='0' cellpadding='5' style='border-collapse:collapse;border-color:#ccc'>"
                 f"<tr style='background:#eee'><th>Method</th><th>{html.escape(res.settings.size_label)}</th>"
                 "<th>Strain e</th><th>R²</th><th>n</th><th>Comments</th></tr>")
        for r in res.methods.values():
            r2 = f"{r.r2:.3f}" if r.r2 is not None and np.isfinite(r.r2) else ""
            com = "; ".join(r.warnings) if r.success else r.message
            grey = " style='color:#999'" if r.extra.get("shape_undetermined") else ""
            H.append(f"<tr{grey}><td>{html.escape(r.name)}</td><td align='right'>{html.escape(r.size_text())}</td>"
                     f"<td align='right'>{html.escape(r.strain_text())}</td><td align='right'>{r2}</td>"
                     f"<td align='right'>{r.n}</td><td style='font-size:11px'>{html.escape(com)}</td></tr>")
        H.append("</table><p style='color:#555;font-size:11px'>Per-reflection methods: mean ± SD between "
                 "reflections. Linear methods: ± standard error of the regression. Upper limits are 2σ. "
                 f"K = {res.settings.K:.4g}; breadths are integral breadths of the instrument-corrected Voigt profile. "
                 "Method assumptions are listed in the Report tab.</p>")
        self.txt_result.setHtml("".join(H))


def main():
    warnings.filterwarnings("ignore", category=UserWarning, module="uncertainties")
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("Völva")
    app.setOrganizationName("UniversityOfSheffield")
    icon = os.path.join(os.path.dirname(os.path.abspath(__file__)), "volva.ico")
    if os.path.isfile(icon):
        app.setWindowIcon(QIcon(icon))
    if sys.platform == "win32":   # own taskbar icon rather than Python's
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("Sheffield.Volva.4")
        except Exception:
            pass
    w = MainWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
