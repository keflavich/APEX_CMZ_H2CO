from astropy import units as u
from astropy import constants
import numpy as np

m_h = constants.m_p

from scipy.optimize import fsolve

def tkin_all(density, sigma, lengthscale, gradient, tdust, crir=1e-17*u.s**-1,
             column=1e22*u.cm**-2, Fx=0*u.erg/u.s/u.cm**2):
    """
    Solve:
        Lamda_(gas-dust) + Lamda_(gas) - Gamma_(turb) - Gamma_(crir) - Gamma_(xray) = 0
    (where any of the heating terms can be set to zero by setting the
    appropriate constant to zero, but the cooling terms cannot be turned off)

    Lambda_(gas-dust) is given by Goldsmith & Langer 1978
    Lambda_(gas) is a parametrization of Goldsmith 2001 by Ao 2013 via Papadapolous 2010 eqn 9
    Gamma_(turb) comes from Pan & Padoan 2009
    Gamma_(crir) is from Tielens 2005
    Gamma_(xray) is direct from Ao 2013; I did not re-examine the derivation of
                 this term

    Uses `scipy.optimize.fsolve` to solve the equation numerically

    Parameters
    ----------
    density : u.cm**-3 equivalent
        Volume density
    sigma : u.km/u.s equivalent
        Line width (sigma, not fwhm)
        (set to zero to turn off turbulent heating)
    lengthscale : u.pc equivalent
        Length scale of turbulent driving
    gradient : u.km/u.s/u.pc equivalent
        Velocity gradient (for line cooling; governs escape fraction)
    tdust : u.K equivalent
        The dust temperature
    crir : u.s**-1 equivalent
        The cosmic ray ionization rate
        (set to zero to turn off cosmic ray heating)
    column : u.cm**-2 equivalent
        The column density of the gas (for X-ray only)
    Fx : u.erg/u.s/u.cm**2 equivalent
        The X-ray flux density
        (set to 0 to turn off X-ray heating)
    """
    assert density.unit.is_equivalent(u.cm**-3)
    assert sigma.unit.is_equivalent(u.km/u.s)
    assert lengthscale.unit.is_equivalent(u.pc)
    assert gradient.unit.is_equivalent(u.km/u.s/u.pc)
    assert crir.unit.is_equivalent(1/u.s)
    N = column.to(u.cm**-2)

    energy_density = u.erg/u.cm**3/u.s


    # c1 = leading constant for gas-dust cooling term
    c1 = 4e-33 * u.cm**3 * u.K**-1.5 *u.erg/u.s
    # c2 = leading constant for gas line cooling term
    c2 = 6e-29 * u.erg/u.cm**3 /u.s *u.K**-3 *(u.km/u.s/u.pc)**-1 * (u.cm**-1.5)**-1
    # ccr = leading term for cosmic ray heating term.  Also, a band.
    ccr = 3.2e-28 * u.erg * u.s**-1

    dvdr = gradient.to(u.km/u.s/u.pc)
    L = lengthscale.to(u.pc)
    n = density.to(u.cm**-3)
    sigma = sigma.to(u.km/u.s)

    # lambda = cooling
    # gamma = heating
    lamgd  = lambda Tk: c1 * n**2 * Tk**(0.5) * (Tk-tdust)
    lamgas = lambda Tk: c2 * n**(1/2.) * Tk**3 * dvdr
    gamturb = n * 2.8 * m_h * (0.5*3**1.5 * sigma**3 / L)
    gamgr = ccr * n * (crir/(1e-17*u.s**-1))
    gammaxray = 1.2e-19 * (n/(1e5*u.cm**-3)) * (Fx/(u.erg*u.cm**-2*u.s**-1)) * (N/(1e22*u.cm**-2))**-0.9 * u.erg/u.s/u.cm**3

    assert gamgr.unit.is_equivalent(energy_density)

    def f(Tk):
        return (lamgd(Tk*u.K).to(u.erg/u.s/u.cm**3).value
                + lamgas(Tk*u.K).to(u.erg/u.s/u.cm**3).value
                - gamgr.value
                - gamturb.to(u.erg/u.s/u.cm**3).value
                - gammaxray.value)

    # 50 is the initial temperature guess
    return fsolve(f, 50)


def tkin_linecooled(density, sigma, lengthscale, gradient, mu=2.35):
    """
    OBSOLETE: use tkin_all
    Equation 21 of Ao et al 2013, with some internal consistency checks.

    Assumes the Goldsmith 2001 cooling rate, i.e.:
    :math:\Lambda_{gas} = 6\\times10^{-29} n^{1/2} T_{kin}^3 \\frac{dv}{dr} [erg s cm^{-3}]:math:

    And uses the Pan and Padoan 2009 heating rate, eqn 10, with an assumed mean heating rate
    :math:\epsilon = 0.5 * \sqrt{3} \sigma_v^3 / L:math:
    (from their Section 1, which gives the sqrt(3) factor that is missing in Ao 2013)

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
    """ Obsolete """
    density = density.to(u.cm**-3).value
    gradient = gradient.to(u.km/u.s/u.pc).value
    lengthscale = lengthscale.to(u.pc).value
    sigma = sigma.to(u.km/u.s).value
    term1 = (16e-8*density**3 + 7920*density**0.5 * gradient * sigma**3 * lengthscale**-1)**0.5 / (12*gradient)
    term2 = 4e-4*density**1.5/(12*gradient)
    return (term1-term2)**(2/3.) * u.K


def analytic_solver():
    """ Obsolete: used to derive analytic solutions when they were possible,
    but that was only for a very limited subset of parameters """
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

