from astropy import units as u
from astropy import constants
import numpy as np

m_h = constants.m_p

def epsilon(x,t):
    '''
    turbulence dissipation rate...
    '''

def turb_heat_rate(density, sigma, mu=2.35, max_length_scale=1*u.pc):
    """
    """
    eps_bar = 0.5 * 3**1.5 * sigma**3 / max_length_scale
    return density * mu * m_h * eps_bar


def tkin_linecooled(density, sigma, lengthscale, gradient, mu=2.35):
    """
    Equation 21 of Ao et al 2013, with some internal consistency checks.

    Assumes the Goldsmith 2001 cooling rate, i.e.:
    :math:\Lambda_{gas} = 6\\times10^{-29} n^{1/2} T_{kin}^3 \\frac{dv}{dr} [erg s cm^{-3}]:math:

    And uses the Pan and Padoan 2009 heating rate, eqn 10, with an assumed mean heating rate
    :math:\epsilon = 0.5 * \sqrt{3} \sigma_v^3 / L:math:

    Parameters
    ----------
    density : u.cm**-3 equivalent
        Volume density
    sigma : u.km/u.s equivalent
        Line width (sigma, not fwhm)

    Example
    -------
    >>> ao_example = tkin_linecooled(10**4.5*u.cm**-3, 20/2.35*u.km/u.s, 5*u.pc, 5*u.km/u.s/u.pc)
    >>> assert np.abs(ao_example-62*u.K) < 1*u.K
    
    """
    assert density.unit.is_equivalent(u.cm**-3)
    assert sigma.unit.is_equivalent(u.km/u.s)
    assert lengthscale.unit.is_equivalent(u.pc)
    assert gradient.unit.is_equivalent(u.km/u.s/u.pc)
    const1 = 6e-29 * u.erg/u.cm**3 /u.s *u.K**-3 *(u.km/u.s/u.pc)**-1 * (u.cm**-1.5)**-1
    const2 = (0.5 * 3**1.5)
    const = const2/const1
    const = (0.5*3**1.5/6e-29) / (u.erg/u.K**3 * u.pc/u.km * u.cm**-1.5)
    assert const2/const1 == const
    turb_heating = const2*density*mu * m_h * sigma**3 / lengthscale
    tkin = (turb_heating / (const1*density**0.5*gradient))**(1/3.)
    return tkin.to(u.K)

def tkin_linedustcooled(density, sigma, lengthscale, gradient):
    density = density.to(u.cm**-3).value
    gradient = gradient.to(u.km/u.s/u.pc).value
    lengthscale = lengthscale.to(u.pc).value
    sigma = sigma.to(u.km/u.s).value
    term1 = (16e-8*density**3 + 7920*density**0.5 * gradient * sigma**3 * lengthscale**-1)**0.5 / (12*gradient)
    term2 = 4e-4*density**1.5/(12*gradient)
    return (term1-term2)**(2/3.) * u.K

from scipy.optimize import fsolve

def tkin_all(density, sigma, lengthscale, gradient, tdust, crir=1e-17*u.s**-1):
    """
    Solve: Lamda_(gas-dust) + Lamda_(gas) - Gamma_(turb) - Gamma_(crir) = 0

    Lambda_(gas-dust) is given by Goldsmith & Langer 1978
    Lambda_(gas) is a parametrization of Goldsmith 2001 by Ao 2013
    Gamma_(turb) comes from Pan & Padoan 2009
    Gamma_(crir) is from Tielens 2005

    Uses `scipy.optimize.fsolve`

    Parameters
    ----------
    density : u.cm**-3 equivalent
        Volume density
    sigma : u.km/u.s equivalent
        Line width (sigma, not fwhm)
    lengthscale : u.pc equivalent
        Length scale of turbulent driving
    gradient : u.km/u.s/u.pc equivalent
        Velocity gradient (for line cooling; governs escape fraction)
    tdust : u.K equivalent
        The dust temperature
    crir : u.s**-1 equivalent
        The cosmic ray ionization rate
    """
    assert density.unit.is_equivalent(u.cm**-3)
    assert sigma.unit.is_equivalent(u.km/u.s)
    assert lengthscale.unit.is_equivalent(u.pc)
    assert gradient.unit.is_equivalent(u.km/u.s/u.pc)
    assert crir.unit.is_equivalent(1/u.s)

    energy_density = u.erg/u.cm**3/u.s


    c1 = 4e-33 * u.cm**3 * u.K**-1.5 *u.erg/u.s
    c2 = 6e-29 * u.erg/u.cm**3 /u.s *u.K**-3 *(u.km/u.s/u.pc)**-1 * (u.cm**-1.5)**-1
    ccr = 3.2e-28 * u.erg * u.s**-1

    dvdr = gradient.to(u.km/u.s/u.pc)
    L = lengthscale.to(u.pc)
    n = density.to(u.cm**-3)
    sigma = sigma.to(u.km/u.s)

    lamgd  = lambda Tk: c1 * n**2 * Tk**(0.5) * (Tk-tdust)
    lamgas = lambda Tk: c2 * n**(1/2.) * Tk**3 * dvdr
    gamturb = n * 2.35 * m_h * (0.5*3**1.5 * sigma**3 / L)
    lamcr = ccr * n * (crir/(1e-17*u.s**-1))

    assert lamcr.unit.is_equivalent(energy_density)

    def f(Tk):
        return (lamgd(Tk*u.K).to(u.erg/u.s/u.cm**3).value
                + lamgas(Tk*u.K).to(u.erg/u.s/u.cm**3).value
                - lamcr.value
                - gamturb.to(u.erg/u.s/u.cm**3).value)

    # 50 is the initial temperature guess
    return fsolve(f, 50)

def solver():
    from sympy import Symbol,solve
    n = Symbol('n')
    Tk = Symbol('Tk')
    Td = Symbol('Td')
    Td = 25
    dvdr = Symbol('dvdr')
    L = Symbol('L')
    c1 = Symbol('c1') # 4e-33
    c2 = Symbol('c2') # 6e-29
    c3 = Symbol('c3') # turb_scalefactor =  (u.cm**-3 * u.g * (u.km/u.s)**3 / u.pc).to(u.erg/u.s/u.cm**3)
    sigma = Symbol('sigma')
    lamgd = c1 * n**2 * Tk**(0.5) * (Tk-25)
    lamgas = c2 * n**(1/2.) * Tk**3 * dvdr
    gamturb = c3 * n * 2.35 * m_h.value * (0.5*3**1.5 * sigma**3 / L)
    gasdustsoln = solve(lamgd+lamgas-gamturb, Tk)
    gassoln = solve(lamgas-gamturb, Tk)
    print "Gas: ",gassoln[0]
    print "Gas+Dust: ",gasdustsoln[0]
    return gassoln, gasdustsoln

