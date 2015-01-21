import pylab as pl
import numpy as np
import aplpy
import wcsaxes
import os
import copy
from astropy.wcs import WCS
from astropy import log
from astropy.io import fits
from astropy import units as u
import paths
from paths import h2copath, figurepath, hpath, fpath
from spectral_cube import SpectralCube, BooleanArrayMask
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits import axisartist as AA

# Close these figures so we can remake them in the appropriate size
for fignum in (1,4,5):
    pl.close(fignum)

recen = dict(x=0.57, y=20e3, width=2.26, height=250e3)
vmin=10
vmax = 200
cmap = pl.cm.RdYlBu_r
figsize = (20,10)

toloop = zip(('TemperatureCube{0}_PiecewiseFromRatio.fits',
              'TemperatureCube_DendrogramObjects{0}.fits',
              'pv_H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e4.fits',
              'pv_H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e5.fits',
              'pv_H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens3e4.fits',
              'pv_H2CO_321220_to_303202{0}_bl_integ_masked_weighted_temperature_dens1e4.fits',
              'pv_H2CO_321220_to_303202{0}_bl_integ_masked_weighted_temperature_dens1e5.fits',
              'pv_H2CO_321220_to_303202{0}_bl_integ_masked_weighted_temperature_dens3e4.fits',),
             ('','_dendro',
              '_directpv1e4','_directpv1e5','_directpv3e4',
              '_directpv1e4_masked','_directpv1e5_masked','_directpv3e4_masked',))

for weight in ("_weight",""):
    for smooth in ("_smooth",""):#"_vsmooth"):
        sncubefile = hpath('APEX_H2CO_303_202{0}_signal_to_noise_cube.fits'.format(smooth))
        sncube = SpectralCube.read(sncubefile)
        snproj = sncube.max(axis=1)
        snproj.wcs.wcs.ctype[0] = 'OFFSET'
        snhdu = snproj.hdu
        snhdu.header['CTYPE1'] = 'GLON'

        fig5 = pl.figure(5, figsize=figsize)
        fig5.clf()
        Fsn = aplpy.FITSFigure(snhdu, convention='calabretta', figure=fig5)
        Fsn.show_grayscale(vmin=0, vmax=10, stretch='linear', invert=True)
        Fsn.add_colorbar()
        Fsn.colorbar.set_axis_label_text('Peak S/N')
        Fsn.recenter(**recen)
        Fsn.tick_labels.set_xformat('d.dd')
        Fsn.save(fpath("big_maps/pv_peaksn{0}.pdf".format(smooth)))

        for ftemplate,outtype in toloop:

            log.info("Starting "+ftemplate.format(smooth))

            if 'pv' not in ftemplate:
                cube = SpectralCube.read(hpath(ftemplate.format(smooth)))
                if weight:
                    wcube = (SpectralCube.read(hpath('APEX_H2CO_303_202_smooth_bl.fits'))
                             if 'smooth' in smooth
                             else SpectralCube.read(hpath('APEX_H2CO_303_202_bl.fits')))
                    mcube = (SpectralCube.read(hpath('APEX_H2CO_303_202_smooth_bl_mask.fits'))
                             if 'smooth' in smooth
                             else SpectralCube.read(hpath('APEX_H2CO_303_202_bl_mask.fits')))
                    if cube.shape != wcube.shape:
                        log.info("Not weighting {0}".format(fn))
                        continue
                    weighted = copy.copy(cube)
                    weighted._data = wcube._data * cube._data
                    mask = (wcube.mask & cube.mask &
                            BooleanArrayMask(mcube._data==1, mcube.wcs))
                    weighted = weighted.with_mask(mask)
                    wcube = wcube.with_mask(mask)

                    pv1 = weighted.sum(axis=1)
                    pv2 = wcube.sum(axis=1)
                    pv = pv1/pv2
                    # Hack to prevent unmatched celestial axis error
                    pv1.wcs.wcs.ctype[0] = 'OFFSET'
                    pv2.wcs.wcs.ctype[0] = 'OFFSET'
                    hdu = copy.copy(pv1.hdu)
                    hdu.data = pv.value
                else:
                    tproj = cube.mean(axis=1)
                    tproj.wcs.wcs.ctype[0] = 'OFFSET'
                    hdu = tproj.hdu
            else:
                hdu = fits.open(hpath(ftemplate.format(smooth)))[0]
            hdu.header['CTYPE1'] = 'GLON'



            # APLPY VERSION
            fig = pl.figure(4, figsize=figsize)
            fig.clf()

            F = aplpy.FITSFigure(hdu,
                                 figure=fig)

            cm = copy.copy(cmap)
            cm.set_bad((0.5,)*3)
            F.show_colorscale(cmap=cm,vmin=vmin,vmax=vmax,
                              aspect=1 if smooth=='' else 2)
            F.recenter(**recen)
            F.tick_labels.set_xformat('d.dd')
            F.add_colorbar()
            F.colorbar.set_axis_label_text('T (K)')
            F.save(fpath("big_maps/pv_tmap{0}{2}{1}.pdf".format(smooth, outtype, weight)))


            color = (0.5,)*3 # should be same as background #888
            F.show_contour(snhdu,
                           levels=[-1,0]+np.logspace(0.20,2).tolist(),
                           colors=([(0.5,0.5,0.5,1)]*2 + 
                                   [color + (alpha,) for alpha in
                                    np.exp(-(np.logspace(0.20,2)-1.7)**2/(2.5**2*2.))]),
                           filled=True,
                           smooth=3,
                           zorder=10, convention='calabretta')


            log.info("Saving " +
                     fpath("big_maps/pv_tmap{0}{2}{1}_masked.pdf".format(smooth,
                                                                         outtype,
                                                                         weight)))
            F.save(fpath("big_maps/pv_tmap{0}{2}{1}_masked.pdf".format(smooth,
                                                                       outtype,
                                                                       weight)))




            # WCSAXES VERSION
            wcs = WCS(hdu.header)

            fig1 = pl.figure(1, figsize=figsize)
            fig1.clf()
            ax = wcsaxes.WCSAxesSubplot(fig1, 1,1,1, wcs=wcs)
            fig1.add_axes(ax)

            ims = ax.imshow(hdu.data, cmap=cm, vmin=vmin, vmax=vmax,
                            aspect=2 if smooth else 1
                           )
            ax.coords[1].set_format_unit(u.km/u.s)

            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            # http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=-0.5 if smooth else 0.05,
                                      axes_class=matplotlib.axes.Axes)

            cb = pl.colorbar(mappable=ims, cax=cax)
            cax.set_ylabel("Temperature [K]")
            cb.ax.yaxis.set_label_position('right')

            ax.set_xlabel("Galactic Longitude")
            ax.set_ylabel("Velocity")
            pl.savefig(fpath("big_maps/pv_wcsaxes_tmap{0}{2}{1}_masked.pdf".format(smooth,
                                                                       outtype,
                                                                       weight)),
                       bbox_inches='tight'
                      )
