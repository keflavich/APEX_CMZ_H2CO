import numpy as np
import pylab as pl
from astropy import units as u
import pyspeckit
from pyspeckit_fitting import (simplemodel, simplemodel2, simple_fitter,
                               simple_fitter2)
from full_cubes import cube_merge_high,pcube_merge_high
from masked_cubes import (cube303, cube303sm, cube303m, cube321m, cube303msm,
                          cube321msm, cube321, cube321sm)
from noise import (noise, noise_cube, sm_noise, cube303nm, cube303nmsm,
                   cube321nm, cube321nmsm)
import paths

def multiscale_fit(pcube, centerx, centery, offset_scale=0.3, savedir=None,
                   savepre=None, exclude=(236,323, 539,654, 972,1049,
                                          1180,1238),
                  guesses=[0.84,37.9,6.5,0.5,0.8,1.0]):

    pcube.specfit.Registry.add_fitter('h2co_simple', simple_fitter2, 6,
                                      multisingle='multi')
    
    pcube.plot_spectrum(centerx, centery)
    pcube.baseline(xtype='pixel', exclude=exclude, order=0, xmin=0,
                   xmax=len(pcube.data), subtract=True,
                   selectregion=True, reset_selection=False,
                   annotate=False)
    pcube.specfit(fittype='h2co_simple', multifit=True,
                      guesses=guesses,
                      limited=[(True,True)] * 6,
                      limits=[(0,20),(-105,225),(1,40),(0,1),(0.3,1.1),(0,1e5)],
                     )
    pcube.baseline(excludefit=True, order=0, subtract=True, annotate=False)
    pcube.specfit(fittype='h2co_simple', multifit=True,
                      guesses=guesses,
                      limited=[(True,True)] * 6,
                      limits=[(0,20),(-105,225),(1,40),(0,1),(0.3,1.1),(0,1e5)],
                     )
    pcube.specfit.plotresiduals(axis=pcube.plotter.axis,
                                      yoffset=-pcube.specfit.parinfo[0].value*offset_scale, clear=False,
                                      color='#444444', label=False)
    pcube.plotter.axis.set_ylim(-3*pcube.specfit.residuals.std()-pcube.specfit.parinfo[0].value*offset_scale,
                                      pcube.plotter.axis.get_ylim()[1])
    pcube.plotter.savefig(paths.fpath('{0}/{1}_r0.pdf'.format(savedir, savepre)))
    r1 = [pcube.specfit.parinfo[3].value]
    er1 = [pcube.specfit.parinfo[3].error]
    sig = [pcube.specfit.parinfo[2].value]
    esig = [pcube.specfit.parinfo[2].error]
    cen = [pcube.specfit.parinfo[1].value]
    ecen = [pcube.specfit.parinfo[1].error]

    for radius in (2,3,4,5,6,7,8,9):
        pcube.plot_apspec([centerx, centery, radius], wunit='pixel')
        pcube.baseline(xtype='pixel', exclude=exclude, order=0, xmin=0,
                             xmax=len(pcube.data), subtract=True,
                             selectregion=True, reset_selection=False,
                             annotate=False)
        pcube.specfit(fittype='h2co_simple', multifit=True,
                          guesses=guesses,
                          limited=[(True,True)] * 6,
                          limits=[(0,20),(-105,225),(1,40),(0,1),(0.3,1.1),(0,1e5)],
                         )
        pcube.baseline(excludefit=True, order=0, subtract=True, annotate=False)
        pcube.specfit(fittype='h2co_simple', multifit=True,
                          guesses=guesses,
                          limited=[(True,True)] * 6,
                          limits=[(0,20),(-105,225),(1,40),(0,1),(0.3,1.1),(0,1e5)],
                         )
        pcube.specfit.plotresiduals(axis=pcube.plotter.axis,
                                          yoffset=-pcube.specfit.parinfo[0].value*offset_scale, clear=False,
                                          color='#444444', label=False)
        pcube.plotter.axis.set_ylim(-3*pcube.specfit.residuals.std()-pcube.specfit.parinfo[0].value*offset_scale,
                                          pcube.plotter.axis.get_ylim()[1])
        pcube.plotter.savefig(paths.fpath('{0}/{1}_r{2}.pdf'.format(savedir,
                                                                    savepre,
                                                                    radius)))
        r1.append(pcube.specfit.parinfo[3].value)
        er1.append(pcube.specfit.parinfo[3].error)
        sig.append(pcube.specfit.parinfo[2].value)
        esig.append(pcube.specfit.parinfo[2].error)
        cen.append(pcube.specfit.parinfo[1].value)
        ecen.append(pcube.specfit.parinfo[1].error)

    pl.figure(2)
    pl.clf()
    ax1 = pl.subplot(2,1,1)
    ax1.errorbar(np.arange(1,10)*7.2, r1, yerr=er1, linestyle='none', marker='s', color='k')
    ax1.set_xlim(0.5*7.2,9.5*7.2)
    #ax1.set_xlabel("Aperture Radius (arcseconds)")
    ax1.set_ylabel("Ratio $3_{2,1}-2_{2,0} / 3_{0,3} - 2_{0,2}$")
    ax1.xaxis.set_ticklabels([])
    ax2 = pl.subplot(2,1,2)
    pl.subplots_adjust(hspace=0)
    ax2.errorbar(np.arange(1,10)*7.2, sig, yerr=esig, linestyle='none', marker='s', color='k')
    ax2.set_xlim(0.5*7.2,9.5*7.2)
    ax2.set_xlabel("Aperture Radius (arcseconds)")
    ax2.set_ylabel("$\sigma$ (km s$^{-1}$)")
    pl.savefig(paths.fpath('{0}/{1}_ratio_vs_scale.pdf'.format(savedir, savepre)))

    return r1,er1,sig,esig,cen,ecen

def multiscale_fit_g08south_hot(center=(434,63)):
    g08south_slice = np.s_[:,center[1]-16:center[1]+16,center[0]-16:center[0]+16]

    g08south_cube = cube_merge_high[g08south_slice]
    g08south_pcube = pyspeckit.Cube(cube=g08south_cube)
    g08south_pcube.xarr.refX = g08south_cube.wcs.wcs.restfrq

    multiscale_fit(g08south_pcube, 16, 16, savedir='g0.8_south', savepre='g0.8_south_hotspot',
                   guesses=[0.5, 40.3, 7.9, 0.57, 0.5, 0.3], offset_scale=0.4)

def multiscale_fit_g08south_cool(center=(426,36)):
    g08south_cool_slice = np.s_[:,center[1]-16:center[1]+16,center[0]-16:center[0]+16]

    g08south_cool_cube = cube_merge_high[g08south_cool_slice]
    g08south_cool_pcube = pyspeckit.Cube(cube=g08south_cool_cube)
    g08south_cool_pcube.xarr.refX = g08south_cool_cube.wcs.wcs.restfrq

    multiscale_fit(g08south_cool_pcube, 16, 16, savedir='g0.8_south', savepre='g0.8_south_coolspot',
                   guesses=[0.5, 50.3, 10.6, 0.3, 0.9, 0.2]+
                           [0.3, 96.3, 3.6, 0.3, 0.9, 0.2]
                   , offset_scale=0.4)

def multiscale_fit_g1pt2_cool(center=(236,95)):
    g1pt2_cool_slice = np.s_[:,center[1]-16:center[1]+16,center[0]-16:center[0]+16]

    g1pt2_cool_cube = cube_merge_high[g1pt2_cool_slice]
    g1pt2_cool_pcube = pyspeckit.Cube(cube=g1pt2_cool_cube)
    g1pt2_cool_pcube.xarr.refX = g1pt2_cool_cube.wcs.wcs.restfrq

    return multiscale_fit(g1pt2_cool_pcube, 16, 16, savedir='g1.2',
                          savepre='g1.2_coolspot', guesses=[0.1, 91, 13.6,
                                                            0.3, 0.9, 0.2],
                          offset_scale=0.4)


def multiscale_fit_g1pt6(center=(53,143)):
    g1pt6_slice = np.s_[:,center[1]-16:center[1]+16,center[0]-16:center[0]+16]

    g1pt6_cube = cube_merge_high[g1pt6_slice]
    g1pt6_pcube = pyspeckit.Cube(cube=g1pt6_cube)
    g1pt6_pcube.xarr.refX = g1pt6_cube.wcs.wcs.restfrq

    return multiscale_fit(g1pt6_pcube, 16, 16, savedir='g1.6',
                          savepre='g1.6_spot',
                          guesses=[0.1, 157, 8, 0.3, 0.9, 0.2]+
                                  [0.1, 60,  9, 0.3, 0.9, 0.2],
                          offset_scale=0.4)

def all_multiscales():
    multiscale_fit_g1pt2_cool()
    multiscale_fit_g08south_hot()
    multiscale_fit_g08south_cool()
    multiscale_fit_g1pt6()
