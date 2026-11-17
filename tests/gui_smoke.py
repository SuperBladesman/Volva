"""Offscreen GUI smoke test: QT_QPA_PLATFORM=offscreen python tests/gui_smoke.py [outdir]"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtWidgets import QApplication, QMessageBox
from PySide6.QtCore import Qt
from volva.gui.app import MainWindow, PlotTab
from volva.io import load_pattern
from volva.pipeline import run_analysis

out = sys.argv[1] if len(sys.argv) > 1 else None
ex = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")
app = QApplication.instance() or QApplication([])
QMessageBox.question = staticmethod(lambda *a, **k: QMessageBox.No)
w = MainWindow(); w.show()
w._load_profile_file(os.path.join(ex, "SYNTHETIC_LaB6_profile.json"))
w.pattern = load_pattern(os.path.join(ex, "SYNTHETIC_CeO2_Lvol34nm_strain0.106pc.xy"))
w._refresh_all()
# run analysis through the worker thread, as the user would
w.analyse()
t0 = time.time()
while w._thread is not None and time.time() - t0 < 120:
    app.processEvents(); time.sleep(0.02)
assert w.result is not None, "analysis did not complete"
print("headline:", w.result.interpretation.headline)
print("table rows:", w.table.rowCount())
# untick first peak and add an hkl label -> re-run
w.table.item(0, 0).setCheckState(Qt.Unchecked)
w.table.item(1, 2).setText("200")
w._rerun()
assert not w.result.peaks[0].use and w.result.peaks[1].hkl == "200"
print("after untick:", w.result.interpretation.headline)
w.copy_methods_text()
assert QApplication.clipboard().text().startswith("Crystallite size and microstrain were determined"), "methods text not copied"
assert "Methods text" in w.txt_result.toPlainText()
# verdict above the headline, tooltips on every tab, keyboard and table stepping
assert w.result.interpretation.verdict and w.result.interpretation.verdict in w.txt_result.toPlainText()
assert all(w.tabs.tabToolTip(i) for i in range(w.tabs.count())), "a tab has no tooltip"
assert "Assumes" in w.tabs.tabToolTip([w.tabs.tabText(i) for i in range(w.tabs.count())].index("Voigt W-H"))
w.peak_index = 0
w._step_peak(1, show=True)
assert w.peak_index == 1 and w.tabs.currentWidget() is w.plot_tabs["Peak"] and w.table.currentRow() == 1
w.table.setCurrentCell(3, 1); app.processEvents()
assert w.peak_index == 3, "moving through the table did not move the Peak tab"
for ax in w.plot_tabs["Pattern"].fig.axes:
    assert "(" not in ax.get_xlabel() or "/" in ax.get_xlabel(), ax.get_xlabel()
for i in range(w.tabs.count()):
    w.tabs.setCurrentIndex(i); app.processEvents()
    if out:
        os.makedirs(out, exist_ok=True)
        w.grab().save(os.path.join(out, f"tab_{i:02d}_{w.tabs.tabText(i).replace(' ', '_')}.png"))
w.cmb_size.setCurrentText("Sphere diameter (K = 4/3)"); w._rerun()
print("sphere:", w.result.interpretation.headline)
print("GUI smoke test passed")
w.close()
