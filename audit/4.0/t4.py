from run_case import *
import patches, synth
def show(tag, D, eps):
    ic, peaks, fits, res, truth = run(D, eps)
    s=summ(res)
    f=lambda v: 'FAIL' if v[0]=='FAIL' else (f"{v[0]:.1f}" + (f"/{v[1]:.3f}%" if v[1] is not None else ""))
    print(f"{tag:22s} D={D} eps={eps}%  " + "  ".join(f"{k[:6]}={f(v)}" for k,v in s.items()), f"| Cag R2={ic.caglioti_r_squared:.3f}")
for D,eps in [(50,0.0),(20,0.0),(100,0.0),(50,0.10)]:
    patches.restore(); show("as-is",D,eps/100)
    patches.fix_fwhm(); show("+FWHM fix",D,eps/100)
    patches.fix_doublet(); show("+FWHM+doublet fix",D,eps/100)
    print()
