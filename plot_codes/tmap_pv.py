import pylab as pl
import numpy as np
import aplpy
import os
import copy
from astropy import log
import paths
from paths import h2copath, figurepath, hpath, fpath
from spectral_cube import SpectralCube
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

# Close these figures so we can remake them in the appropriate size
for fignum in (4,5):
    pl.close(fignum)

cmap = pl.cm.RdYlBu_r
figsize = (20,10)

for smooth in ("","_smooth",):#"_vsmooth"):
    sncubefile = hpath('APEX_H2CO_303_202{0}_signal_to_noise_cube.fits'.format(smooth))
    sncube = SpectralCube.read(sncubefile)
    snproj = sncube.max(axis=1)
    snproj.wcs.wcs.ctype[0] = 'OFFSET'
    snhdu = snproj.hdu
    snhdu.header['CTYPE1'] = 'GLON'

    for ftemplate,outtype in zip(('TemperatureCube_PiecewiseFromRatio.fits',
                                  'TemperatureCube_DendrogramObjects{0}_Piecewise.fits'),
                                 ('','dendro')):

        log.info("Starting "+ftemplate.format(smooth))
        cube = SpectralCube.read(hpath(ftemplate.format(smooth)))
        tproj = cube.mean(axis=1)
        tproj.wcs.wcs.ctype[0] = 'OFFSET'
        hdu = tproj.hdu
        hdu.header['CTYPE1'] = 'GLON'

        fig5 = pl.figure(5, figsize=figsize)
        fig5.clf()
        Fsn = aplpy.FITSFigure(snhdu, convention='calabretta', figure=fig5)
        Fsn.show_grayscale(vmin=0, vmax=10, stretch='linear', invert=True)
        Fsn.add_colorbar()
        Fsn.colorbar.set_axis_label_text('Peak S/N')
        Fsn.recenter(0.55, 50e3, width=2.3, height=300e3)
        Fsn.tick_labels.set_xformat('d.dd')
        Fsn.save(fpath("big_maps/pv_peaksn{0}_{1}.pdf".format(smooth, outtype)))


        fig = pl.figure(4, figsize=figsize)
        fig.clf()

        F = aplpy.FITSFigure(hdu,
                             convention='calabretta',
                             figure=fig)

        cm = copy.copy(cmap)
        cm.set_bad((0.5,)*3)
        F.show_colorscale(cmap=cm,vmin=15,vmax=200)
        F.recenter(0.55, 50e3, width=2.3, height=300e3)
        F.tick_labels.set_xformat('d.dd')
        F.add_colorbar()
        F.colorbar.set_axis_label_text('T (K)')
        F.save(fpath("big_maps/pv_tmap{0}_{1}.pdf".format(smooth, outtype)))


        color = (0.5,)*3 # should be same as background #888
        F.show_contour(snhdu,
                       levels=[-1,0]+np.logspace(0.20,2).tolist(),
                       colors=([(0.5,0.5,0.5,1)]*2 + 
                               [color + (alpha,) for alpha in
                                np.exp(-(np.logspace(0.20,2)-1.7)**2/(2.5**2*2.))]),
                       filled=True,
                       smooth=3,
                       zorder=10, convention='calabretta')


        log.info("Saving "+fpath("big_maps/pv_tmap{0}_{1}_masked.pdf".format(smooth, outtype)))
        F.save(fpath("big_maps/pv_tmap{0}_{1}_masked.pdf".format(smooth, outtype)))
