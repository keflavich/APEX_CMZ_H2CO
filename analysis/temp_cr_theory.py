import numpy as np
from astropy import units as u

def cr_cooling(density, gradient, crcoolingrate, dustdom=True):
    if dustdom:
        return ((1.6e-7 * density**3 + 768. * density**0.5*gradient *
                 (crcoolingrate/(1e-17/u.s)))**0.5 / (12*gradient) -
                (4e-4*density**1.5 / (12*gradient))**(2./3.))
    else: # line dominated
        return ((6*(density/(10**4.5/u.cm**3))**(1/6.) *
                 (gradient/(5*u.km/u.s/u.pc))**(-1/3.) *
                 (crcoolingrate/(1e-17/u.s))**(1/3.)))

def gas_cooling(density, gastem, dusttem):
    return 4e-33 * density**2 * gastem**0.5 * (gastem-dusttem)

