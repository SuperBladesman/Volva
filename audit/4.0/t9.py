from run_case import *
import patches, synth
orig=synth.inst_profile
patches.fix_fwhm(); patches.fix_doublet()
for label,fn in [("instrument eta~0.5 (as before)",orig),
                 ("purely Gaussian instrument", lambda tt:(orig(tt)[0]*1.3, 1e-6))]:
    synth.inst_profile=fn
    row=[]
    for D in [20,50,100,200]:
        ic,peaks,fits,res,truth=run(D,0); row.append(f"{D}->{res['Scherrer']['average_size_nm']:.0f}")
    print(f"{label:32s} Scherrer (FWHM+doublet bugs fixed): "+", ".join(row))
