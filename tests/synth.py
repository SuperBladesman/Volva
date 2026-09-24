"""Synthetic Bragg-Brentano patterns with known size/strain, true Voigt convolution, Ka doublet, Poisson noise."""
import numpy as np
from scipy.special import voigt_profile
from itertools import product
LAM1, LAM2, R21 = 0.1540562, 0.1544390, 0.5
K_TRUE = 0.94

def voigt_fwhm(fg, fl):  # Olivero & Longbothum
    return 0.5346*fl + np.sqrt(0.2166*fl**2 + fg**2)

def refl_cubic(a_nm, rule, tt_max=125):
    out = {}
    for h,k,l in product(range(0,9), repeat=3):
        if (h,k,l)==(0,0,0) or not (h>=k>=l): continue
        if not rule(h,k,l): continue
        d = a_nm/np.sqrt(h*h+k*k+l*l)
        s = LAM1/(2*d)
        if s>=1: continue
        tt = 2*np.degrees(np.arcsin(s))
        if tt>tt_max: continue
        key=round(tt,3)
        perms = len(set(__import__('itertools').permutations((h,k,l))))
        signs = 2**sum(1 for v in (h,k,l) if v)
        out[key] = out.get(key,0) + perms*signs
    return sorted(out.items())

def inst_profile(tt):
    """Instrument FWHM (deg) components: Gaussian via Caglioti, Lorentzian small linear-in-tan term."""
    th = np.radians(tt/2)
    fg = np.sqrt(0.0030*np.tan(th)**2 - 0.0015*np.tan(th) + 0.0025)
    fl = 0.020 + 0.015*np.tan(th)
    return fg, fl

def pattern(refls, D_nm=None, eps=0.0, tt=np.arange(10,125,0.02), scale=2e4, bg=150, seed=0, lam2=True, asym=0.0):
    rng = np.random.default_rng(seed)
    y = np.full_like(tt, bg, dtype=float) + 30*np.exp(-(tt-10)/30)
    truth = []
    for pos, mult in refls:
        th = np.radians(pos/2)
        fg_i, fl_i = inst_profile(pos)
        Dv = D_nm(pos) if callable(D_nm) else D_nm
        fl_s = np.degrees(K_TRUE*LAM1/(Dv*np.cos(th))) if Dv else 0.0   # Lorentzian size
        fg_s = np.degrees(4*eps*np.tan(th))                                   # Gaussian strain
        fg = np.sqrt(fg_i**2 + fg_s**2); fl = fl_i + fl_s
        amp = scale*mult/(np.sin(th)**2*np.cos(th))*np.exp(-2*(np.sin(th)/LAM1)**2*0.01)/50
        sig = fg/(2*np.sqrt(2*np.log(2))); gam = fl/2
        comps = [(pos, 1.0, LAM1)]
        if lam2:
            p2 = np.degrees(2*np.arcsin(np.sin(th)*LAM2/LAM1)); comps.append((p2, R21, LAM2))
        a = asym(pos) if callable(asym) else asym      # axial-divergence-like low-angle tail
        for p, r, _ in comps:
            sc = np.where(tt < p, 1.0 + a, 1.0 - a)     # split profile: area, height and FWHM unchanged
            y += amp*r*voigt_profile((tt-p)/sc, sig, gam)*0.05
        truth.append(dict(tt=pos, fwhm_tot=voigt_fwhm(fg,fl), fwhm_inst=voigt_fwhm(fg_i,fl_i),
                          fl_size=fl_s, fg_strain=fg_s))
    y = rng.poisson(np.maximum(y,0)).astype(float)
    return tt, y, truth

LAB6 = lambda: refl_cubic(0.415689, lambda h,k,l: True)
CEO2 = lambda: refl_cubic(0.54113, lambda h,k,l: (h%2==k%2==l%2))


def write_xy(path, tt, y, header=None):
    with open(path, "w") as f:
        if header:
            f.write(f"# {header}\n")
        for a, b in zip(tt, y):
            f.write(f"{a:.4f} {b:.1f}\n")
