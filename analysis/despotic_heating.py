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
import despotic
import os
import numpy as np

# Use the Milky Way GMC file as a base
gmc=cloud('cloud.desp')

# from despotic.chemistry import NL99
# gmc.setChemEq(network=NL99)

def tkin_all(density, sigma, lengthscale, gradient, tdust, crir=1e-17*u.s**-1,
             ISRF=1, tdust_rad=None, turbulence=True):

    assert density.unit.is_equivalent(u.cm**-3)
    assert sigma.unit.is_equivalent(u.km/u.s)
    assert lengthscale.unit.is_equivalent(u.pc)
    assert gradient.unit.is_equivalent(u.km/u.s/u.pc)
    assert crir.unit.is_equivalent(1/u.s)

    gmc.sigmaNT = sigma.to(u.cm/u.s).value
    gmc.Td = tdust.to(u.K).value
    gmc.rad.TradDust = gmc.Td if tdust_rad is None else tdust_rad.to(u.K).value
    gmc.dVdr = gradient.to(u.s**-1).value
    gmc.rad.chi = ISRF

    # These are both per hydrogen, but we want to specify per particle, and
    # we're assuming the particles are H2
    gmc.rad.ionRate = crir.to(u.s**-1).value * 2
    gmc.nH = density.to(u.cm**-3).value * 2

    def turb_heating(cloud, driving_scale=lengthscale):
        """ Turbulent heating rate depends on cloud linewidth
        (sigma_nonthermal) and driving scale of the turbulence """
        if turbulence:
            gamturb = 1.4 * constants.m_p * cloud.nH*u.cm**-3 * (0.5*3**1.5 * (cloud.sigmaNT*u.cm/u.s)**3 / (driving_scale))
            return [(gamturb/(cloud.nH*u.cm**-3)).to(u.erg/u.s).value, 0]
        else:
            return [0,0]

    gmc.setTempEq(escapeProbGeom='LVG', PsiUser=turb_heating)
    #energy_balance = gmc.dEdt()

    return gmc.Tg

if __name__ == "__main__":
    import matplotlib
    matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

    densities = np.logspace(3,7,20)
    tem = [tkin_all(n*u.cm**-3, 10*u.km/u.s, lengthscale=5*u.pc,
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
    linewidths = np.logspace(np.log10(0.5), np.log10(30), 15)
    tem2 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem3 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem4 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                    crir=1e-14*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem5 = [tkin_all(1e6*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem10 = [tkin_all(1e6*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                    crir=1e-14*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem6 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=1*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem7 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                    crir=1e-14*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem8 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                     gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=25*u.K,
                     crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem9 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                     tdust_rad=10*u.K,
                     gradient=20*u.km/u.s/u.pc, tdust=25*u.K,
                     crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem11 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                     tdust_rad=10*u.K,
                     gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     crir=2e-14*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem12 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                     tdust_rad=10*u.K,
                     gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     crir=1e-15*u.s**-1) for sigma in ProgressBar(linewidths)]

    FWHM = np.sqrt(8*np.log(2))
    fig = pl.figure(2)
    pl.clf()
    ax = pl.gca()
    ax.plot(linewidths*FWHM, tem2,  'k--', alpha=0.5, linewidth=2,
            label='$\zeta_{CR}=1e-17$ s$^{-1}$\n $n=10^4$ cm$^{-3}$\n'
                  '$L=5$ pc\n $dv/dr=5$ km/s/pc\n'
                  '$T_D=25$K\n $T_D(rad)=10$K')
    ax.plot(linewidths*FWHM, tem7,  'k:',  alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$, $n=10^4$')
    ax.plot(linewidths*FWHM, tem9,  'k-',  alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-17$, $n=10^4$ $dv/dr=20$')
    ax.plot(linewidths*FWHM, tem3,  'r--', alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-17$, $n=10^5$')
    ax.plot(linewidths*FWHM, tem6,  'r-',  alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-17$, $n=10^5$ L=1 pc')
    ax.plot(linewidths*FWHM, tem4,  'r:',  alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$, $n=10^5$')
    ax.plot(linewidths*FWHM, tem10, 'b:',  alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$, $n=10^6$')
    ax.plot(linewidths*FWHM, tem5,  'b--', alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-17$, $n=10^6$')
    ax.plot(linewidths*FWHM, tem8,  'r-.', alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-17$, $n=10^5$ $T_D(rad)=25$ K')
    ax.plot(linewidths*FWHM, tem11, 'r-', alpha=0.5, linewidth=4, label='$\zeta_{CR}=1e-13$, $n=10^5$')
    ax.plot(linewidths*FWHM, tem12, 'r-', alpha=0.5, linewidth=3, label='$\zeta_{CR}=1e-15$, $n=10^5$')
    ax.set_xlabel("Line FWHM (km s$^{-1}$)")
    ax.set_ylabel("Temperature (K)")
    ax.set_ylim(0,150)
    ax.set_xlim(0,linewidths.max()*FWHM)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    ax.legend(loc='center left', fontsize=16, bbox_to_anchor=(1.0, 0.55))
    pl.savefig(paths.fpath("despotic/TvsSigma.png"), bbox_inches='tight')
    pl.savefig(paths.fpath("despotic/TvsSigma.pdf"), bbox_inches='tight')

    from astropy.table import Table
    pcfittable = Table.read(paths.apath('fitted_line_parameters_Chi2Constraints.ipac'),
                            format='ascii.ipac')

    lolim = pcfittable['tmax1sig_chi2'] > 340
    maps = np.char.startswith(pcfittable['Source_Name'], 'Map')
    ok = ~np.isnan(pcfittable['tmin1sig_chi2']) & (pcfittable['width'] < 40) & (pcfittable['h2coratio321303']/pcfittable['eh2coratio321303'] > 5) & pcfittable['is_good'].astype('bool')
    flags = {'is_map': maps,
             'is_lolim': lolim,
             'is_ok': ok}
    # Don't plot these for now...
    pcfittable = pcfittable[(~lolim) & ok]
    maps = np.char.startswith(pcfittable['Source_Name'], 'Map')
    lolim_conservative = pcfittable['tmax1sig_chi2'] > 150

    mask = maps&~lolim_conservative
    ax.errorbar(pcfittable['width'][mask]*(8*np.log(2))**0.5,
                 pcfittable['temperature_chi2'][mask],
                 yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[mask],
                       (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[mask]],
                 capsize=0,
                 markersize=10,
                 markeredgecolor='none',
                 linestyle='none', marker='s', linewidth=0.5, alpha=0.6, color='r')
    mask = maps&lolim_conservative
    ax.plot(pcfittable['width'][mask]*(8*np.log(2))**0.5,
             pcfittable['tmin1sig_chi2'][mask],
             marker='^',
             markersize=10,
             markeredgecolor='none',
             color='r',
             alpha=0.4,
             linestyle='none')

    ax.set_xlabel("Line FWHM (km s$^{-1}$)")
    ax.set_ylabel("Temperature (K)")
    ax.set_ylim(0,150)
    fig.savefig(paths.fpath('despotic/chi2_temperature_vs_linewidth_byfield.pdf'),
                             bbox_inches='tight')

    mask = (~maps)&(~lolim_conservative)
    ax.errorbar(pcfittable['width'][mask]*(8*np.log(2))**0.5,
                 pcfittable['temperature_chi2'][mask],
                 yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[mask],
                       (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[mask]],
                 capsize=0,
                 markeredgecolor='none',
                 markersize=10,
                 linestyle='none', marker='s', linewidth=0.5, alpha=0.6, color='b')

    mask = (~maps)&lolim_conservative
    ax.plot(pcfittable['width'][mask]*(8*np.log(2))**0.5,
             pcfittable['tmin1sig_chi2'][mask],
             marker='^',
             markersize=10,
             markeredgecolor='none',
             color='b',
             alpha=0.4,
             linestyle='none')

    ax.set_ylim(0,150)
    fig.savefig(paths.fpath('despotic/chi2_temperature_vs_linewidth_fieldsandsources.pdf'),
                             bbox_inches='tight')


    pl.draw(); pl.show()
