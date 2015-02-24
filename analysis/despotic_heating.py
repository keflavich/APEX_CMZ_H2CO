"""
Copied from despotic/examples/gmcChem and modified
"""
import matplotlib
import pylab as pl
from astropy import units as u
from astropy import constants
import paths
from astropy.utils.console import ProgressBar

# Import the despotic library and the NL99 network; also import numpy
from despotic import cloud
#from despotic.chemistry import NL99
import despotic
import os
import numpy as np

# Use the Milky Way GMC file as a base
gmc=cloud('cloud.desp')

def tkin_all(density, sigma, lengthscale, gradient, tdust, crir=1e-17*u.s**-1,
             ISRF=1):

    assert density.unit.is_equivalent(u.cm**-3)
    assert sigma.unit.is_equivalent(u.km/u.s)
    assert lengthscale.unit.is_equivalent(u.pc)
    assert gradient.unit.is_equivalent(u.km/u.s/u.pc)
    assert crir.unit.is_equivalent(1/u.s)

    gmc.sigmaNT = sigma.to(u.cm/u.s).value
    gmc.Td = gmc.rad.TradDust = tdust.to(u.K).value
    gmc.dVdr = gradient.to(u.s**-1).value
    gmc.rad.chi = ISRF

    # These are both per hydrogen, but we want to specify per particle, and
    # we're assuming the particles are H2
    gmc.rad.ionRate = crir.to(u.s**-1).value * 2
    gmc.nH = density.to(u.cm**-3).value * 2

    def turb_heating(cloud, driving_scale=lengthscale):
        """ Turbulent heating rate depends on cloud linewidth (sigma_nonthermal) and driving scale of the turbulence """
        gamturb = 1.4 * constants.m_p * cloud.nH*u.cm**-3 * (0.5*3**1.5 * (cloud.sigmaNT*u.cm/u.s)**3 / (driving_scale))
        return [(gamturb/(cloud.nH*u.cm**-3)).to(u.erg/u.s).value, 0]

    gmc.setTempEq(escapeProbGeom='LVG', PsiUser=turb_heating)
    #energy_balance = gmc.dEdt()

    return gmc.Tg

if __name__ == "__main__":

    densities = np.logspace(3,7,20)
    tem = [tkin_all(n*u.cm**-3, 10*u.km/u.s, lengthscale=10*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                    crir=1e-17*u.s**-1) for n in ProgressBar(densities)]
    pl.figure(1)
    pl.clf()
    pl.plot(densities, tem, 'k--', label='CRIR=1e-17, $\sigma=10$ km/s')
    pl.xlabel(r'$\log\,N_{\rm H}$')
    pl.ylabel('Temperature (K)')
    pl.xscale('log')
    pl.legend(loc='best')
    pl.savefig(paths.fpath("despotic/TvsN.png"))


    linewidths = np.arange(0.5,30,2)
    tem2 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=10*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem3 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=10*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem4 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=10*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                    crir=1e-14*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem5 = [tkin_all(1e6*u.cm**-3, sigma*u.km/u.s, lengthscale=10*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem6 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=1*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem7 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=10*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                    crir=1e-14*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem8 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=10*u.pc,
                     gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     ISRF=1000,
                     crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem9 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=10*u.pc,
                    gradient=20*u.km/u.s/u.pc, tdust=25*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]

    FWHM = np.sqrt(8*np.log(2))
    pl.figure(2)
    pl.clf()
    pl.plot(linewidths*FWHM, tem2, 'k--', alpha=0.5, linewidth=2, label='CRIR=1e-17, $n=10^4$')
    pl.plot(linewidths*FWHM, tem7, 'k:',  alpha=0.5, linewidth=2, label='CRIR=1e-14, $n=10^4$')
    pl.plot(linewidths*FWHM, tem9, 'k-',  alpha=0.5, linewidth=2, label='CRIR=1e-17, $n=10^4$ $dv/dr=20$')
    pl.plot(linewidths*FWHM, tem3, 'r--', alpha=0.5, linewidth=2, label='CRIR=1e-17, $n=10^5$')
    pl.plot(linewidths*FWHM, tem6, 'r-',  alpha=0.5, linewidth=2, label='CRIR=1e-17, $n=10^5$ L=1 pc')
    pl.plot(linewidths*FWHM, tem4, 'r:',  alpha=0.5, linewidth=2, label='CRIR=1e-14, $n=10^5$')
    pl.plot(linewidths*FWHM, tem5, 'b:',  alpha=0.5, linewidth=2, label='CRIR=1e-14, $n=10^6$')
    pl.plot(linewidths*FWHM, tem8, 'r-.', alpha=0.5, linewidth=2, label='CRIR=1e-17, $n=10^5$ ISRF=1000')
    pl.xlabel(r'FWHM $= 2.35 \sigma$ km s$^{-1}$')
    pl.ylabel('Temperature (K)')
    pl.ylim(0,150)
    pl.legend(loc='best')
    pl.savefig(paths.fpath("despotic/TvsSigma.png"))

    pl.draw(); pl.show()
