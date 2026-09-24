from run_case import *
import synth
ic, peaks, fits, res, truth = run(50, 0.0)
print("STANDARD: fitted vs true instrument FWHM")
for p in ic.standard_peaks:
    th=np.radians(p['two_theta']/2); fg,fl=synth.inst_profile(p['two_theta']); tr=synth.voigt_fwhm(fg,fl)
    f=p['fit_result']
    print(f"{p['two_theta']:8.3f} fit={p['fwhm']:.4f} true={tr:.4f} err={100*(p['fwhm']/tr-1):+5.1f}%  eta={f.fraction:.2f} doublet={f.doublet_fitted} 2sigma={2*f.sigma:.4f} ({100*(2*f.sigma/tr-1):+.1f}%)")
print("SAMPLE:")
tmap={round(t['tt'],1):t for t in truth}
for p in peaks:
    t=min(truth,key=lambda t:abs(t['tt']-p.two_theta_deg))
    f=min([f for f in fits if f.success],key=lambda f:abs(f.center_deg-p.two_theta_deg))
    th=p.theta_rad
    D_ideal=0.94*synth.LAM1/(np.radians(t['fl_size'])*np.cos(th))
    D_v=0.94*synth.LAM1/(p.beta_corrected_rad*np.cos(th))
    print(f"{p.two_theta_deg:8.3f} obs={p.fwhm_obs_deg:.4f} true={t['fwhm_tot']:.4f} ({100*(p.fwhm_obs_deg/t['fwhm_tot']-1):+5.1f}%) 2sig={2*f.sigma:.4f} eta={f.fraction:.2f} inst={p.fwhm_inst_deg:.4f}/true {t['fwhm_inst']:.4f}  D_volva={D_v:6.1f} D_true={D_ideal:.1f}")
