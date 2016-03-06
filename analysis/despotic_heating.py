"""
Copied from despotic/examples/gmcChem and modified
"""
from __future__ import print_function
import matplotlib
import pylab as pl
pl.switch_backend('Qt4Agg')
from astropy import units as u
from astropy import constants
from astropy.table import Table
import paths
from paths import fpath
from astropy.utils.console import ProgressBar
import pprint

# Import the despotic library and the NL99 network; also import numpy
from despotic import cloud
import despotic
import os
import numpy as np

# Use the Milky Way GMC file as a base
gmc=cloud('cloud.desp')

from despotic.chemistry import NL99
# gmc.setChemEq(network=NL99)

def turb_heating_generator(lengthscale=1*u.pc, turbulence=True):
    def turb_heating(cloud, lengthscale=lengthscale):
        """ Turbulent heating rate depends on cloud linewidth
        (sigma_nonthermal) and driving scale of the turbulence
        DESPOTIC wants units of erg/s/H (per hydrogen), so the turbulent
        heating rate n sigma^3 / L is divided by n to get just sigma^3/L
        """
        if turbulence:
            gamturb = (1.4 * constants.m_p *
                       (0.5*3**1.5 * (cloud.sigmaNT*u.cm/u.s)**3 / (lengthscale)))
            return [(gamturb).to(u.erg/u.s).value, 0]
        else:
            return [0,0]
    return turb_heating

def fiducial_case(sigma=5.0*u.km/u.s, tdust=25*u.K, tdust_rad=10*u.K,
                  gradient=5*u.km/u.s/u.pc, ISRF=0, crir=1e-14*u.s**-1,
                  turbulence=False, density=10**4*u.cm**-3,
                  lengthscale=1*u.pc):
    print()
    print("Fiducial case from section 5.1 para 4: ")
    cloud
    gmc.sigmaNT = sigma.to(u.cm/u.s).value
    gmc.Td = tdust.to(u.K).value
    gmc.rad.TradDust = gmc.Td if tdust_rad is None else tdust_rad.to(u.K).value
    gmc.dVdr = gradient.to(u.s**-1).value
    gmc.rad.chi = ISRF

    # These are both per hydrogen, but we want to specify per particle, and
    # we're assuming the particles are H2
    gmc.rad.ionRate = crir.to(u.s**-1).value * 2
    gmc.nH = density.to(u.cm**-3).value * 2

    turb_heating = turb_heating_generator(lengthscale, turbulence=turbulence)

    print("Before: ",)
    pprint.pprint(gmc.dEdt(), width=1)
    print(gmc.setTempEq(escapeProbGeom='LVG', PsiUser=turb_heating))
    print("After: ",)
    pprint.pprint(gmc.dEdt(), width=1)

    return gmc.Tg


def tkin_all(density, sigma, lengthscale, gradient, tdust, crir=1e-17*u.s**-1,
             ISRF=1, tdust_rad=None, turbulence=True, gmc=gmc, reload_gmc=True,
             chemistry=False):

    assert density.unit.is_equivalent(u.cm**-3)
    assert sigma.unit.is_equivalent(u.km/u.s)
    assert lengthscale.unit.is_equivalent(u.pc)
    assert gradient.unit.is_equivalent(u.km/u.s/u.pc)
    assert crir.unit.is_equivalent(1/u.s)

    if reload_gmc:
        gmc=cloud('cloud.desp')

    gmc.sigmaNT = sigma.to(u.cm/u.s).value
    gmc.Td = tdust.to(u.K).value
    gmc.rad.TradDust = gmc.Td if tdust_rad is None else tdust_rad.to(u.K).value
    gmc.dVdr = gradient.to(u.s**-1).value
    gmc.rad.chi = ISRF

    # These are both per hydrogen, but we want to specify per particle, and
    # we're assuming the particles are H2
    gmc.rad.ionRate = crir.to(u.s**-1).value * 2
    gmc.nH = density.to(u.cm**-3).value * 2

    turb_heating = turb_heating_generator(lengthscale, turbulence=turbulence)

    try:
        gmc.setTempEq(escapeProbGeom='LVG', PsiUser=turb_heating)
    except despotic.despoticError as ex:
        print(ex)
        return np.nan

    if chemistry:

        gmc.setChemEq(network=NL99)
        gmc.setTempEq(escapeProbGeom='LVG', PsiUser=turb_heating)

    return gmc.Tg

def case_study(row, gmc=gmc):
    beam_pc = (30/(np.sqrt(8*np.log(2)))*u.arcsec*8.5*u.kpc).to(u.pc,
                                                                u.dimensionless_angles())
    print("Case study for object ID ",row['_idx'])
    gf=row['gausscorrfactor']
    print("Gaussian correction factor: ",gf)
    print("Density: ",10**row['density_chi2']*u.cm**-3)
    print("Line width: ",row['v_rms']*u.km/u.s*gf)
    print("Lengthscale: ",2*((row['reff']*u.pc*gf)**2-(beam_pc)**2)**0.5)
    print("Tdust: ",row['higaldusttem']*u.K)
    print("Tdust,rad: ",(row['higaldusttem']*u.K *
                         (1-np.exp(-(10**row['logh2column']/1e24)))))

    lengthscale = 2*((row['reff']*u.pc*gf)**2-(beam_pc)**2)**0.5

    T0 = tkin_all(density=10**row['density_chi2']*u.cm**-3,
                      sigma=row['v_rms']*u.km/u.s*gf,
                      lengthscale=lengthscale,
                      gradient=5*u.km/u.s/u.pc, #min(5,row['v_rms']/row['reff'])*u.km/u.s/u.pc,
                      tdust=row['higaldusttem']*u.K,
                      crir=1e-17*u.s**-1,
                      ISRF=1,
                      tdust_rad=(row['higaldusttem']*u.K *
                                 (1-np.exp(-(10**row['logh2column']/1e24)))),
                  gmc=gmc)
    print("Initial temperature: ",T0)
    cool0 = gmc.dEdt(PsiUser=turb_heating_generator(lengthscale))
    print("CO cooling: ",cool0['LambdaLine']['co'])
    print("O cooling: ",cool0['LambdaLine']['o'])
    print("C cooling: ",cool0['LambdaLine']['c'])
    print("C+ cooling: ",cool0['LambdaLine']['c+'])
    print("oH2 cooling: ",cool0['LambdaLine']['oh2'])
    print("pH2 cooling: ",cool0['LambdaLine']['ph2'])
    print("HD cooling: ",cool0['LambdaLine']['hd'])
    gmc.setChemEq(network=NL99)

    T1 = tkin_all(density=10**row['density_chi2']*u.cm**-3,
                      sigma=row['v_rms']*u.km/u.s*gf,
                      lengthscale=2*((row['reff']*u.pc*gf)**2-(beam_pc)**2)**0.5,
                      gradient=5*u.km/u.s/u.pc, #min(5,row['v_rms']/row['reff'])*u.km/u.s/u.pc,
                      tdust=row['higaldusttem']*u.K,
                      crir=1e-17*u.s**-1,
                      ISRF=1,
                      tdust_rad=(row['higaldusttem']*u.K *
                                 (1-np.exp(-(10**row['logh2column']/1e24)))),
                  gmc=gmc)
    print("Chemical equilbrium temperature: ",T1)
    cool1 = gmc.dEdt(PsiUser=turb_heating_generator(lengthscale))
    print("CO cooling: ",cool1['LambdaLine']['co'])
    print("O cooling: ",cool1['LambdaLine']['o'])
    print("C cooling: ",cool1['LambdaLine']['c'])
    print("C+ cooling: ",cool1['LambdaLine']['c+'])
    print("oH2 cooling: ",cool1['LambdaLine']['oh2'])
    print("pH2 cooling: ",cool1['LambdaLine']['ph2'])
    print("HD cooling: ",cool1['LambdaLine']['hd'])

    print()
    print("The same, but with an enhanced IRSF = 1000x local")
    T0 = tkin_all(density=10**row['density_chi2']*u.cm**-3,
                      sigma=row['v_rms']*u.km/u.s*gf,
                      lengthscale=2*((row['reff']*u.pc*gf)**2-(beam_pc)**2)**0.5,
                      gradient=5*u.km/u.s/u.pc, #min(5,row['v_rms']/row['reff'])*u.km/u.s/u.pc,
                      tdust=row['higaldusttem']*u.K,
                      crir=1e-17*u.s**-1,
                      ISRF=1000,
                      tdust_rad=(row['higaldusttem']*u.K *
                                 (1-np.exp(-(10**row['logh2column']/1e24)))),
                  gmc=gmc)
    print("Initial temperature: ",T0)
    cool0 = gmc.dEdt(PsiUser=turb_heating_generator(lengthscale))
    print("CO cooling: ",cool0['LambdaLine']['co'])
    print("O cooling: ",cool0['LambdaLine']['o'])
    print("C cooling: ",cool0['LambdaLine']['c'])
    print("C+ cooling: ",cool0['LambdaLine']['c+'])
    print("oH2 cooling: ",cool0['LambdaLine']['oh2'])
    print("pH2 cooling: ",cool0['LambdaLine']['ph2'])
    print("HD cooling: ",cool0['LambdaLine']['hd'])
    gmc.setChemEq(network=NL99)

    T1 = tkin_all(density=10**row['density_chi2']*u.cm**-3,
                      sigma=row['v_rms']*u.km/u.s*gf,
                      lengthscale=2*((row['reff']*u.pc*gf)**2-(beam_pc)**2)**0.5,
                      gradient=5*u.km/u.s/u.pc, #min(5,row['v_rms']/row['reff'])*u.km/u.s/u.pc,
                      tdust=row['higaldusttem']*u.K,
                      crir=1e-17*u.s**-1,
                      ISRF=1000,
                      tdust_rad=(row['higaldusttem']*u.K *
                                 (1-np.exp(-(10**row['logh2column']/1e24)))),
                  gmc=gmc)
    print("Chemical equilbrium temperature: ",T1)
    cool1 = gmc.dEdt(PsiUser=turb_heating_generator(lengthscale))
    print("CO cooling: ",cool1['LambdaLine']['co'])
    print("O cooling: ",cool1['LambdaLine']['o'])
    print("C cooling: ",cool1['LambdaLine']['c'])
    print("C+ cooling: ",cool1['LambdaLine']['c+'])
    print("oH2 cooling: ",cool1['LambdaLine']['oh2'])
    print("pH2 cooling: ",cool1['LambdaLine']['ph2'])
    print("HD cooling: ",cool1['LambdaLine']['hd'])


if __name__ == "__main__":
    import matplotlib
    matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

    densities = np.logspace(3,7,20)

    linewidths = np.arange(0.5,30,2)
    linewidths = np.logspace(np.log10(0.5), np.log10(30), 15)
    # do tem4 first because it crashes sometimes
    tem4 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                     gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                     crir=1e-14*u.s**-1,
                     dampfactor=0.02) for sigma in ProgressBar(linewidths)]

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


    fiducial = tem2 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem3 = [tkin_all(1e5*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                    gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     tdust_rad=10*u.K,
                    crir=1e-17*u.s**-1) for sigma in ProgressBar(linewidths)]
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
    tem6[-2] = np.nan # this point is bad
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
    tem11 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                     tdust_rad=10*u.K,
                     gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     crir=2e-14*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem12 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                     tdust_rad=10*u.K,
                     gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     crir=1e-15*u.s**-1) for sigma in ProgressBar(linewidths)]
    tem14 = [tkin_all(1e4*u.cm**-3, sigma*u.km/u.s, lengthscale=5*u.pc,
                     tdust_rad=10*u.K,
                     gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                     crir=1e-14*u.s**-1,
                     turbulence=False) for sigma in ProgressBar(linewidths)]


    linewidths2 = np.linspace(1,linewidths.max(),50.)
    def rho(sig, rho_5kms=10**4.25*u.cm**-3):
        return rho_5kms *(sig/(5*u.km/u.s))**-2
    def L(sig, L_5kms=5*u.pc):
        return L_5kms *(sig/(5*u.km/u.s))**0.7
    tem13 = []
    tem15 = []
    lambdaline = []
    lambdadust = []
    gammaturb = []
    for sigma in ProgressBar(linewidths2):
        tem13.append(tkin_all(density=rho(sigma*u.km/u.s),
                      sigma=sigma*u.km/u.s,
                      lengthscale=L(sigma*u.km/u.s),
                      tdust_rad=10*u.K,
                      gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                      crir=1e-17*u.s**-1))
        dedt = gmc.dEdt(PsiUser=turb_heating_generator(L(sigma*u.km/u.s)))
        lambdaline.append(np.sum(dedt['LambdaLine'].values()))
        gammaturb.append(dedt['PsiUserGas'])
        lambdadust.append(dedt['LambdaDust'])
        tem15.append(tkin_all(density=1e4*u.cm**-3,
                      sigma=sigma*u.km/u.s,
                      lengthscale=L(sigma*u.km/u.s),
                      tdust_rad=10*u.K,
                      gradient=5*u.km/u.s/u.pc, tdust=25*u.K,
                      crir=1e-17*u.s**-1))


    print("Polynomial approximations for the fiducial case: ")
    print("Linear (bad fit): T = sigma*{0} + {1}".format(*np.polyfit(linewidths, fiducial, 1)))
    print("2nd order (good fit): T = sigma**2*{0} + sigma*{1} + {2}".format(*np.polyfit(linewidths, fiducial, 2)))

    # Plotting starts here
    FWHM = np.sqrt(8*np.log(2))
    fig = pl.figure(2)
    pl.clf()
    ax = pl.gca()
    ax.plot(linewidths*FWHM, fiducial,  'k--', alpha=0.5, linewidth=2,
            label='$\zeta_{CR}=1e-17$ s$^{-1}$\n $n=10^4$ cm$^{-3}$\n'
                  '$L=5$ pc\n $dv/dr=5$ km/s/pc\n'
                  '$T_D=25$K\n $T_D(rad)=10$K')
    ax.plot(linewidths*FWHM, tem7,  'k:',  zorder=-5, alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$')
    ax.plot(linewidths*FWHM, tem4,  'r:',  zorder=-5, alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$, $n=10^5$')
    ax.plot(linewidths*FWHM, tem10, 'b:',  zorder=-5, alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$, $n=10^6$')
    ax.plot(linewidths*FWHM, tem11, 'k-',  zorder=-5, alpha=0.2, linewidth=6, label='$\zeta_{CR}=2e-14$')
    ax.plot(linewidths*FWHM, tem12, 'k-',  zorder=-5, alpha=0.3, linewidth=4, label='$\zeta_{CR}=1e-15$')
    ax.plot(linewidths*FWHM, tem9,  'k-',  zorder=-5, alpha=0.5, linewidth=2, label='$dv/dr=20$')
    ax.plot(linewidths*FWHM, tem3,  'r--', zorder=-5, alpha=0.5, linewidth=2, label='$n=10^5$')
    ax.plot(linewidths*FWHM, tem5,  'b--', zorder=-5, alpha=0.5, linewidth=2, label='$n=10^6$')
    ax.plot(linewidths*FWHM, tem6,  'r-',  zorder=-5, alpha=0.5, linewidth=2, label='$n=10^5$ $L=1$')
    ax.plot(linewidths*FWHM, tem8,  'r-.', zorder=-5, alpha=0.5, linewidth=2, label='$n=10^5$ $T_D(rad)=25$')
    ax.plot(linewidths2*FWHM, tem13, 'g-', zorder=-5, alpha=0.5, linewidth=2, label=r'$n=10^{4.25}\sigma_5^{-2}$, $L=5\sigma_5^{0.7}$')
    ax.plot(linewidths2*FWHM, tem15, 'g--', zorder=-5, alpha=0.5, linewidth=2, label=r'$L=5\sigma_5^{0.7}$')
    ax.plot(linewidths*FWHM, tem14, 'g:',  zorder=-5, alpha=0.5, linewidth=2, label=r'$\zeta_{CR}=1e-14$, no turbulence')
    ax.set_xlabel("Line FWHM (km s$^{-1}$)")
    ax.set_ylabel("Temperature (K)")
    ax.set_ylim(0,150)
    ax.set_xlim(2,linewidths.max()*FWHM)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    ax.legend(loc='center left', fontsize=16, bbox_to_anchor=(1.0, 0.55))
    pl.savefig(paths.fpath("despotic/TvsSigma.png"), bbox_inches='tight')
    pl.savefig(paths.fpath("despotic/TvsSigma.pdf"), bbox_inches='tight')

    def plot_one_component(linewidths, data, markerline, alpha=0.5,
                           linewidth=2, label='', figname=None,
                           include_fiducial=False, **kwargs):
        fig3 = pl.figure(3)
        fig3.clf()
        ax3 = fig3.gca()
        ax3.plot(linewidths, data, markerline,  zorder=-5, alpha=alpha, linewidth=linewidth, label=label, **kwargs)
        if include_fiducial:
            ax3.plot(linewidths, fiducial, 'k--',  zorder=-10, alpha=0.5, linewidth=2, label='Fiducial')
        ax3.set_xlabel("Line FWHM (km s$^{-1}$)")
        ax3.set_ylabel("Temperature (K)")
        ax3.set_ylim(0,150)
        ax3.set_xlim(2,linewidths.max())
        if figname is not None:
            fig3.savefig(figname, bbox_inches='tight')

    plot_one_component(linewidths*FWHM, fiducial,  'k--', alpha=0.5, linewidth=2,
                       label='$\zeta_{CR}=1e-17$ s$^{-1}$\n $n=10^4$ cm$^{-3}$\n'
                             '$L=5$ pc\n $dv/dr=5$ km/s/pc\n'
                             '$T_D=25$K\n $T_D(rad)=10$K',
                      figname=paths.fpath("despotic/TvsSigma_fiducial.pdf"))
    plot_one_component(linewidths*FWHM, tem14, 'g:',  alpha=0.5, linewidth=2, label=r'$\zeta_{CR}=1e-14$, no turbulence',
                       figname=paths.fpath("despotic/TvsSigma_CRonly.pdf"))
    plot_one_component(linewidths2*FWHM, tem15, 'g--', alpha=0.5, linewidth=2, label=r'$L=5\sigma_5^{0.7}$',
                       figname=paths.fpath("despotic/TvsSigma_sizelinewidthonly.pdf"))
    plot_one_component(linewidths*FWHM, tem7,  'k:',  alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$',
                       figname=paths.fpath("despotic/TvsSigma_CRm14_only.pdf"))
    plot_one_component(linewidths*FWHM, tem4,  'r:',  alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$, $n=10^5$',
                       figname=paths.fpath("despotic/TvsSigma_CRm14_n5_only.pdf"))
    plot_one_component(linewidths*FWHM, tem10, 'b:',  alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$, $n=10^6$',
                       figname=paths.fpath("despotic/TvsSigma_CRm14_n6_only.pdf"))
    plot_one_component(linewidths*FWHM, tem11, 'k-',  alpha=0.2, linewidth=6, label='$\zeta_{CR}=2e-14$',
                       figname=paths.fpath("despotic/TvsSigma_CR2m14_only.pdf"))
    plot_one_component(linewidths*FWHM, tem12, 'k-',  alpha=0.3, linewidth=4, label='$\zeta_{CR}=1e-15$',
                       figname=paths.fpath("despotic/TvsSigma_CRm15_only.pdf"))
    plot_one_component(linewidths*FWHM, tem9,  'k-',  alpha=0.5, linewidth=2, label='$dv/dr=20$',
                       figname=paths.fpath("despotic/TvsSigma_dvdr20_only.pdf"))
    plot_one_component(linewidths*FWHM, tem3,  'r--', alpha=0.5, linewidth=2, label='$n=10^5$',
                       figname=paths.fpath("despotic/TvsSigma_n5_only.pdf"))
    plot_one_component(linewidths*FWHM, tem5,  'b--', alpha=0.5, linewidth=2, label='$n=10^6$',
                       figname=paths.fpath("despotic/TvsSigma_n6_only.pdf"))
    plot_one_component(linewidths*FWHM, tem6,  'r-',  alpha=0.5, linewidth=2, label='$n=10^5$ $L=1$',
                       figname=paths.fpath("despotic/TvsSigma_n5_L1_only.pdf"))
    plot_one_component(linewidths*FWHM, tem8,  'r-.', alpha=0.5, linewidth=2, label='$n=10^5$ $T_D(rad)=25$',
                       figname=paths.fpath("despotic/TvsSigma_n5_td25_only.pdf"))
    plot_one_component(linewidths2*FWHM, tem13, 'g-', alpha=0.5, linewidth=2, label=r'$n=10^{4.25}\sigma_5^{-2}$, $L=5\sigma_5^{0.7}$',
                       figname=paths.fpath("despotic/TvsSigma_isobar_L0.7_only.pdf"))


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




    from dendrograms import (catalog, catalog_sm, dend, dendsm)
    smooth=''
    cat = catalog
    sn = (cat['ratio321303']/cat['eratio321303'])
    sngt50 = sn > 50
    sn25_50 = (sn > 25) & (sn < 50)
    ok = (np.isfinite(sn) & (cat['Stot321'] < cat['Stot303']) & ~(cat['bad']) &
          (cat['Smean321'] > 0) &
          (cat['e321'] > 0) &
          (~cat['IsNotH2CO']) & (~cat['IsAbsorption']))
    ok = np.array(ok, dtype='bool')
    gt5 = (sn>5)

    hot = cat['temperature_chi2'] > 150
    #gcorfactor = gaussian_correction.gaussian_correction(catalog['Smin303']/catalog['Smax303'])
    gcorfactor = cat['gausscorrfactor']
    masks = (gt5 & ~sngt50 & ~sn25_50 & ok,
             sn25_50 & gt5 & ok,
             sngt50 & gt5 & ok,
             ok & ~gt5)
    is_leaf = np.array(cat['is_leaf'])
    leaf_masks = [np.array(mm, dtype='bool') for mask in masks for mm in (mask & is_leaf, mask & ~is_leaf)]
    # mask1 & leaf, mask1 & not leaf, mask2 & leaf, mask2 & not leaf....
    # Make the not-leaves be half as bright
    masks_colors = zip(leaf_masks,
                       ('b','b','g','g','r','r',    'k','k'),
                       (0.5,0.2, 0.6,0.3, 0.7,0.35, 0.3,0.15),
                       (8,7,9,8,10,9,5,4),
                      )
    pl.figure(12).clf()
    fig12, ax12 = pl.subplots(num=12)
    ax12.errorbar(cat['v_rms'][hot]*np.sqrt(8*np.log(2))*gcorfactor[hot], [149]*hot.sum(),
                  lolims=True, linestyle='none', capsize=0, alpha=0.3,
                  marker='^', color='r')
    for mask,color,alpha,markersize in masks_colors:
        ax12.errorbar(cat['v_rms'][mask]*np.sqrt(8*np.log(2))*gcorfactor[mask], cat['temperature_chi2'][mask],
                      #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                      markersize=10 if any(mask & is_leaf) else 5,
                      markeredgecolor='none',
                      linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax12.set_xlabel(r"Line FWHM (km s$^{-1}$)")
        ax12.set_ylabel("Temperature (K)")

    ax12.plot(linewidths*FWHM, fiducial,  'k--', alpha=0.5, linewidth=2,
              label='$\zeta_{CR}=1e-17$ s$^{-1}$\n $n=10^4$ cm$^{-3}$\n'
                    '$L=5$ pc\n $dv/dr=5$ km/s/pc\n'
                    '$T_D=25$K\n $T_D(rad)=10$K')
    ax12.plot(linewidths*FWHM, tem7,  'k:',  zorder=-5, alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$')
    ax12.plot(linewidths*FWHM, tem9,  'k-',  zorder=-5, alpha=0.5, linewidth=2, label='$dv/dr=20$')
    ax12.plot(linewidths*FWHM, tem11, 'k-',  zorder=-5, alpha=0.2, linewidth=6, label='$\zeta_{CR}=2e-14$')
    ax12.plot(linewidths*FWHM, tem12, 'k-',  zorder=-5, alpha=0.3, linewidth=4, label='$\zeta_{CR}=1e-15$')
    ax12.plot(linewidths*FWHM, tem3,  'r--', zorder=-5, alpha=0.5, linewidth=2, label='$n=10^5$')
    ax12.plot(linewidths*FWHM, tem6,  'r-',  zorder=-5, alpha=0.5, linewidth=2, label='$n=10^5$ $L=1$')
    ax12.plot(linewidths*FWHM, tem4,  'r:',  zorder=-5, alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$, $n=10^5$')
    ax12.plot(linewidths*FWHM, tem10, 'b:',  zorder=-5, alpha=0.5, linewidth=2, label='$\zeta_{CR}=1e-14$, $n=10^6$')
    ax12.plot(linewidths*FWHM, tem5,  'b--', zorder=-5, alpha=0.5, linewidth=2, label='$n=10^6$')
    ax12.plot(linewidths*FWHM, tem8,  'r-.', zorder=-5, alpha=0.5, linewidth=2, label='$n=10^5$ $T_D(rad)=25$')
    ax12.plot(linewidths2*FWHM, tem13, 'g-', zorder=-5, alpha=0.5, linewidth=2, label=r'$n=10^{4.25}\sigma_5^{-2}$, $L=5\sigma_5^{0.7}$')
    ax12.plot(linewidths2*FWHM, tem15, 'g--', zorder=-5, alpha=0.5, linewidth=2, label=r'$L=5\sigma_5^{0.7}$')
    ax12.plot(linewidths*FWHM, tem14, 'g:',  zorder=-5, alpha=0.5, linewidth=2, label=r'$\zeta_{CR}=1e-14$, no turbulence')

    ax12.set_xlim([2,70])
    ax12.set_ylim([0,150])
    fig12.savefig(fpath('despotic/temperature_vs_rmsvelocity{0}.pdf'.format(smooth)))
    wide = cat['v_rms']*gcorfactor > 48/np.sqrt(8*np.log(2))
    ax12.errorbar([50.5] * (wide & is_leaf).sum(),
                  cat['temperature_chi2'][wide&is_leaf],
                  lolims=True, linestyle='none', capsize=0, alpha=0.3,
                  markersize=10,
                  marker='>', color='r')
    ax12.errorbar([50.5] * (wide & ~is_leaf).sum(),
                  cat['temperature_chi2'][wide&(~is_leaf)],
                  lolims=True, linestyle='none', capsize=0, alpha=0.1,
                  markersize=5,
                  marker='>', color='r')
    ax12.set_xlim([2,25])
    fig12.savefig(fpath('despotic/temperature_vs_rmsvelocity_xzoom{0}.pdf'.format(smooth)))


    brick_sw_id = ((catalog['x_cen'] - 0.241)**2 + (catalog['y_cen']+0.0057)**2).argmin()
    print("Brick SW: {0}".format(brick_sw_id))
    print(catalog['_idx','Smean303','ratio321303','higaldusttem','tmin1sig_chi2','temperature_chi2','tmax1sig_chi2'][brick_sw_id])
    row = catalog[brick_sw_id]
    case_study(row)

    print()
    print()

    brick_ne_id = ((catalog['x_cen'] - 0.2615)**2 + (catalog['y_cen']+0.0283)**2).argmin()
    print("Brick NE: {0}".format(brick_ne_id))
    print(catalog['_idx','Smean303','ratio321303','higaldusttem','tmin1sig_chi2','temperature_chi2','tmax1sig_chi2'][brick_ne_id])
    row = catalog[brick_ne_id]
    case_study(row)

    fiducial_case()

    pl.draw(); pl.show()
