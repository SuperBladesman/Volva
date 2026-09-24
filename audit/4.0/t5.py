from run_case import *
import patches, synth, voigt_decon as vd
patches.fix_fwhm(); patches.fix_doublet()
lam=synth.LAM1
for D,eps in [(20,0),(50,0),(100,0),(200,0),(50,0.001)]:
    ic, peaks, fits, res, truth = run(D, eps)
    # instrument G, L vs 2theta from standard fits
    tt=np.array([p['two_theta'] for p in ic.standard_peaks]); th=np.radians(tt/2)
    GL=np.array([vd.decompose(p['fwhm'], p['fit_result'].fraction) for p in ic.standard_peaks])
    A=np.column_stack([np.tan(th)**2,np.tan(th),np.ones_like(th)]); cG=np.linalg.lstsq(A,GL[:,0]**2,rcond=None)[0]
    B=np.column_stack([np.tan(th),1/np.cos(th)]); cL=np.linalg.lstsq(B,GL[:,1],rcond=None)[0]
    Ds=[];xs=[];ys=[]
    for p in peaks:
        f=min([f for f in fits if f.success],key=lambda f:abs(f.center_deg-p.two_theta_deg))
        Go,Lo=vd.decompose(p.fwhm_obs_deg,f.fraction); t=p.theta_rad
        Gi=np.sqrt(max(cG@[np.tan(t)**2,np.tan(t),1],0)); Li=cL@[np.tan(t),1/np.cos(t)]
        Ls=max(Lo-Li,1e-6); Gs=np.sqrt(max(Go**2-Gi**2,0))
        Ds.append(0.94*lam/(np.radians(Ls)*np.cos(t)))
        xs.append(4*np.sin(t)); ys.append(np.radians(Gs)*np.cos(t))
    eps_fit = np.linalg.lstsq(np.array(xs)[:,None], np.array(ys), rcond=None)[0][0]
    s=summ(res)
    print(f"D_true={D:4d} eps_true={eps*100:.2f}% | Voigt-decon size (Lorentzian part) = {np.mean(Ds):6.1f} ± {np.std(Ds):4.1f} nm, strain(Gaussian part)={eps_fit*100:.3f}%  | Volva Scherrer={s['Scherrer'][0]:.1f}, H-W={s['Halder-Wagner'][0]:.1f}/{s['Halder-Wagner'][1]:.3f}%")
