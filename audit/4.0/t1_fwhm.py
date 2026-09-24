# Test 1: does _pseudo_voigt_fwhm match the true FWHM of lmfit's PseudoVoigtModel?
import numpy as np, volva_engine as ve
from lmfit.models import PseudoVoigtModel
m=PseudoVoigtModel()
x=np.linspace(-5,5,200001)
print("lmfit fwhm expr:", m.make_params()['fwhm'].expr)
for eta in [0,0.25,0.5,0.75,1.0]:
    y=m.eval(x=x,amplitude=1,center=0,sigma=0.1,fraction=eta)
    true=ve.PeakFitter._numerical_fwhm(x,y)
    volva=ve.PeakFitter._pseudo_voigt_fwhm(0.1,eta)
    print(f"eta={eta:.2f} true FWHM={true:.5f}  Volva={volva:.5f}  error={100*(volva/true-1):+.1f}%")
