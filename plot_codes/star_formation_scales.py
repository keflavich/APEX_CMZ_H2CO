import numpy as np
import paths
from astropy import units as u
from astropy import constants
from astropy import constants as c

def sonicmass(temperature=10*u.K,
              vrms=None,
              cloud_mach=10.,
              particle_mass=2.8*constants.m_p,
              surfacedensity=1e21*u.cm**-2, toomreQ=1,
              scale_height=10*u.pc,
              p=2,
              eos_gamma=1.):
    """
    Parameters
    ----------
    eos_gamma : float
        The equation of state parameter gamma.  1 = isothermal, 5/3 = adiabatic
    """
    c_s = (eos_gamma * constants.k_B * temperature / particle_mass)**0.5

    if vrms is not None:
        cloud_mach = (vrms / c_s * np.sqrt(3)).decompose()

    #m = 2**(1.5) * toomreQ * c_s**4 / constants.G**2 / (surfacedensity * particle_mass)
    R_sonic = scale_height * cloud_mach**(-2/(p-1.))
    m = 2/3. * c_s**2 * R_sonic / constants.G
    return m.to(u.M_sun)

def jeansmass(temperature=10*u.K, muh2=2.8, density=1e4*u.cm**-3):
    m = (((2 * c.k_B * temperature)**1.5 / ((muh2*c.m_p*c.G)**1.5 *
                                            (muh2*density*c.m_p)**(0.5)))).cgs
    return m.to(u.M_sun)


def approx_temperature(sigma):
    """ 2nd order poly approximation from despotic_heating using the fiducial
    relation """
    sigma = sigma.to(u.km/u.s).value
    return (sigma**2*0.102957453893 + sigma*4.535474454 + 4.65792876797)*u.K

if __name__ == "__main__":
    import pylab as pl

    linewidths = np.linspace(1.0/2.35,100/2.35,500)
    sonic_masses = sonicmass(temperature=approx_temperature(linewidths*u.km/u.s), vrms=linewidths*u.km/u.s)
    sonic_masses_high = sonicmass(temperature=approx_temperature(linewidths*u.km/u.s), vrms=linewidths*u.km/u.s, scale_height=20*u.pc)
    sonic_masses_pfivethirds = sonicmass(temperature=approx_temperature(linewidths*u.km/u.s), vrms=linewidths*u.km/u.s, p=5/3.)
    jeans_masses = jeansmass(temperature=approx_temperature(linewidths*u.km/u.s))
    jeans_masses_high = jeansmass(temperature=approx_temperature(linewidths*u.km/u.s), density=1e5*u.cm**-3)

    fig = pl.figure(4)
    pl.clf()
    ax = pl.gca()
    ax.loglog(linewidths*2.35, sonic_masses, linewidth=3, alpha=0.5, label="$M_{sonic}$")
    ax.loglog(linewidths*2.35, sonic_masses_high, linewidth=3, alpha=0.5, label="$M_{sonic}$ $h=20$")
    ax.loglog(linewidths*2.35, sonic_masses_pfivethirds, linewidth=3, alpha=0.5, label="$M_{sonic}$ $p=5/3$")
    ax.loglog(linewidths*2.35, jeans_masses, linestyle='--', linewidth=3, alpha=0.5, label="$M_{Jeans,thermal}$")
    ax.loglog(linewidths*2.35, jeans_masses_high, linestyle='--', linewidth=3, alpha=0.5, label="$M_{Jeans,thermal}$ $n=10^5$")
    ax.set_ylim(0.1,100)
    
    ax.set_xlabel("Line FWHM (km s$^{-1}$)")
    ax.set_ylabel("Mass ($M_\odot$)")
    ax.legend(loc='best')
    fig.savefig(paths.fpath("despotic/sonic_jeans_vs_linewidth.pdf"), bbox_inches='tight')
