import numpy as np, volva_engine as ve
_orig_fwhm = ve.PeakFitter._pseudo_voigt_fwhm
_orig_fit_peak = ve.PeakFitter.fit_peak
def fix_fwhm():
    ve.PeakFitter._pseudo_voigt_fwhm = staticmethod(lambda sigma, fraction: 2.0*sigma)
def fix_doublet():
    def fit_peak(self, two_theta, intensity, peak_guess_deg, fit_window_deg=2.5, enable_doublet=True):
        half=fit_window_deg/2; m=(two_theta>=peak_guess_deg-half)&(two_theta<=peak_guess_deg+half)
        x=two_theta[m]; y=intensity[m]
        if len(x)<ve.MIN_FIT_POINTS: return ve.PeakFitResult(success=False,message='few pts')
        amp=np.max(y)
        # width guess from data half-max crossings
        above=x[y>=amp/2]; w=max(above[-1]-above[0], 2*np.median(np.diff(x))) if len(above)>1 else 0.1
        sg=w/2
        d=self._delta_2theta_doublet(peak_guess_deg) if enable_doublet else 0
        try:
            if enable_doublet: return self._fit_doublet(x,y,peak_guess_deg,amp,sg,d)
            return self._fit_singlet(x,y,peak_guess_deg,amp,sg)
        except Exception as e:
            return ve.PeakFitResult(success=False,message=str(e))
    ve.PeakFitter.fit_peak = fit_peak
def restore():
    ve.PeakFitter._pseudo_voigt_fwhm = staticmethod(_orig_fwhm); ve.PeakFitter.fit_peak=_orig_fit_peak
