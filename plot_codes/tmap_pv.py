import pylab as pl
import numpy as np
import aplpy
import os
import copy
from astropy import log
from astropy.io import fits
import paths
from paths import h2copath, figurepath, hpath, fpath
from spectral_cube import SpectralCube
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

# Close these figures so we can remake them in the appropriate size
for fignum in (4,5):
    pl.close(fignum)

recen = dict(x=0.57, y=20e3, width=2.26, height=250e3)
vmax = 150
cmap = pl.cm.RdYlBu_r
figsize = (20,10)

for smooth in ("","_smooth",):#"_vsmooth"):
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

    for ftemplate,outtype in zip(('TemperatureCube{0}_PiecewiseFromRatio.fits',
                                  'TemperatureCube_DendrogramObjects{0}.fits',
                                  'pv_H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e4.fits',
                                  'pv_H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e5.fits',
                                  'pv_H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens3e4.fits',
                                 ),
                                 ('','_dendro','_directpv1e4','_directpv1e5','_directpv3e4')):

        log.info("Starting "+ftemplate.format(smooth))

        if 'pv' not in ftemplate:
            cube = SpectralCube.read(hpath(ftemplate.format(smooth)))
            tproj = cube.mean(axis=1)
            tproj.wcs.wcs.ctype[0] = 'OFFSET'
            hdu = tproj.hdu
        else:
            hdu = fits.open(hpath(ftemplate.format(smooth)))[0]
        hdu.header['CTYPE1'] = 'GLON'



        fig = pl.figure(4, figsize=figsize)
        fig.clf()

        F = aplpy.FITSFigure(hdu,
                             convention='calabretta',
                             figure=fig)

        cm = copy.copy(cmap)
        cm.set_bad((0.5,)*3)
        F.show_colorscale(cmap=cm,vmin=15,vmax=vmax)
        F.recenter(**recen)
        F.tick_labels.set_xformat('d.dd')
        F.add_colorbar()
        F.colorbar.set_axis_label_text('T (K)')
        F.save(fpath("big_maps/pv_tmap{0}{1}.pdf".format(smooth, outtype)))


        color = (0.5,)*3 # should be same as background #888
        F.show_contour(snhdu,
                       levels=[-1,0]+np.logspace(0.20,2).tolist(),
                       colors=([(0.5,0.5,0.5,1)]*2 + 
                               [color + (alpha,) for alpha in
                                np.exp(-(np.logspace(0.20,2)-1.7)**2/(2.5**2*2.))]),
                       filled=True,
                       smooth=3,
                       zorder=10, convention='calabretta')


        log.info("Saving "+fpath("big_maps/pv_tmap{0}{1}_masked.pdf".format(smooth, outtype)))
        F.save(fpath("big_maps/pv_tmap{0}{1}_masked.pdf".format(smooth, outtype)))
