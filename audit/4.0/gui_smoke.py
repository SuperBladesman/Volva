import tkinter as tk, traceback, numpy as np
import volva_gui as g
root=tk.Tk(); app=g.VolvaApp(root)
# simulate loading a sample without file dialog
import volva_engine as ve
t,i,m=ve.XRDLoader.load('sample_xrd_data.txt')
g.filedialog.askopenfilename=lambda **k: 'sample_xrd_data.txt'
try:
    app._load_sample(); print("load sample OK")
except Exception as e:
    print("load sample FAILED:", type(e).__name__, e)
try:
    app.sample_two_theta=t; app.sample_intensity=i
    app._run_peak_processing(); print("peaks:", len(app.analysed_peaks))
    app._run_all(); print("run all OK:", list(app.analysis_results))
    for v in ['xrd_bg','xrd_fits','gallery','inspector','scherrer_d','wh','hw','ssp','monshi','comparison']:
        app._switch_plot(v)
    print("all views rendered")
except Exception:
    traceback.print_exc()
root.destroy()
