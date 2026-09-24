import numpy as np, volva_engine as ve, synth, logging, warnings
warnings.filterwarnings('ignore'); logging.disable(logging.WARNING)
def run(D, eps, snip=40, doublet=True, step=0.02, prom=0.05, seed=1, verbose=False, fw=2.5):
    tt = np.arange(10,125,step)
    wl = ve.WAVELENGTHS['Cu Ka']
    fitter = ve.PeakFitter(wl)
    t_s, y_s, _ = synth.pattern(synth.LAB6(), None, 0, tt=tt, seed=seed)
    ic = ve.InstrumentalCorrection()
    msg = ic.process_standard(t_s, y_s, fitter, enable_doublet=doublet, min_r_squared=0.7, prominence_factor=prom)
    t, y, truth = synth.pattern(synth.CEO2(), D, eps, tt=tt, seed=seed+1)
    peaks, fits, log = ve.PeakProcessor.process_sample(t, y, fitter, inst_correction=ic, wavelength_nm=wl['ka1'],
                         prominence_factor=prom, enable_doublet=doublet, snip_iterations=snip, fit_window_deg=fw)
    res = ve.run_all_methods(peaks, wl['ka1'], 0.94)
    if verbose: print(msg); print(log)
    return ic, peaks, fits, res, truth
def summ(res):
    out={}
    for k,r in res.items():
        if not r.get('success'): out[k]=('FAIL',None); continue
        out[k]=(r.get('crystallite_size_nm', r.get('average_size_nm')), r.get('microstrain_percent'))
    return out
