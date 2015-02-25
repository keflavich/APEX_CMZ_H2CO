"""
Python implementation of Rosolowsky & Blitz 2006 Appendix B
"""
import numpy as np
from numpy import exp,pi,sqrt
from scipy.special import erf as errorf

def gaussian_correction(clip_fraction, npow=4, reference_clip=0.05):
    """
    Parameters
    ----------
    clip_fraction : float
        The ratio Tmin / Tmax, i.e. the fraction of peak at which you have
        clipped.  Must be <1
    npow : int
        Should be 4 for variance...
    reference_clip : float
        The reference clipping value (as a fraction of peak) to calibrate to
    """

    # Rosolowsky '06 uses xref=2, which corresponds to clipping at 0.135335
    xref = np.sqrt(2*np.log(1./reference_clip))
    xmax = np.sqrt(2*np.log(1./clip_fraction))

    numerator = gint_eval(xref, npow)/gint_eval(xref, npow-2)
    denominator = gint_eval(xmax, npow)/gint_eval(xmax, npow-2)
    correction = numerator/denominator
    return correction


def gint_eval(x, npow):

    if npow == 0:
      return errorf(x/sqrt(2.))
    elif npow == 1:
      return 1-exp(-x**2./2.)
    elif npow == 2:
      return errorf(x/sqrt(2.))-sqrt(2./pi)*x*exp(-x**2./2.)
    elif npow == 3:
      return 1-exp(-x**2./2.)*(x**2./2.+1)
    elif npow == 4:
      return errorf(x/sqrt(2.))-(x+x**3./3.)*sqrt(2./pi)*exp(-x**2./2.)
    elif npow == 5:
      return 1-exp(-x**2./2.)*(1+x**2./2.+x**4./8.)
    elif npow == 6:
      return errorf(x/sqrt(2.))-sqrt(2./pi)*exp(-x**2./2.)*(x+x**3./3.+x**5./15.)
    else:
        raise ValueError("npow must be 0..5")

