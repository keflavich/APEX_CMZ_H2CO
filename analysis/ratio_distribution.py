from scipy.integrate import quad
import numpy as np

def ratio_distribution(z, mux, muy, sigmax, sigmay):
    """ z is the probability? """
    def a(z):
        return (z**2/sigmax**2 + 1./sigmay**2)**0.5
    def b(z):
        return mux/sigmax**2 * z + muy/sigmay**2
    def c(z):
        return np.exp(0.5*b(z)**2/a(z)**2 - 0.5*(mux**2/sigmax**2 + muy**2/sigmay**2))
    def phi(z):
        def integrand(u):
            return (2*np.pi)**-0.5 * np.exp(-0.5*u**2)
        return quad(integrand, -np.infty, z)[0]

    term1 = b(z) * c(z) / a(z)**3 * (2*np.pi)**-0.5 / sigmax / sigmay * (2*phi(b(z)/a(z)) - 1)
    term2 = (a(z)**2 * np.pi * sigmax * sigmay)**-1 * np.exp(-0.5 * ((mux/sigmax)**2 + (muy/sigmay)**2))

    return term1 + term2

def ratio_error(mux, muy, sigmax, sigmay):
    """
    Approximation that is "good if y unlikely to be negative" (it can't be in
    our distributions, so that's good)
    """
    pass
