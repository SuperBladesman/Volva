from run_case import *
import synth
for D in [5,10,20]:
  for step in [0.02,0.01]:
    for snip in [40,150]:
        ic, peaks, fits, res, truth = run(D, 0, step=step, snip=snip)
        s=summ(res)
        # background error: compare SNIP bg to true bg at a peak
        tt=np.arange(10,125,step); _,y,_=synth.pattern(synth.CEO2(),D,0,tt=tt,seed=2)
        sub,bg=ve.BackgroundProcessor.subtract(tt,y,n_iterations=snip)
        truebg=150+30*np.exp(-(tt-10)/30)
        i=np.argmin(abs(tt-28.55))
        print(f"D={D:3d} step={step} SNIP={snip:3d}: peaks={len(peaks):2d} Scherrer={s['Scherrer'][0] if s['Scherrer'][0]!='FAIL' else 'FAIL'}  H-W={s['Halder-Wagner'][0]}  | bg at (111) peak: SNIP={bg[i]:.0f} true={truebg[i]:.0f}, peak height={y[i]-truebg[i]:.0f}")
