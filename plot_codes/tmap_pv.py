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

for ftemplate,outtype in zip(('TemperatureCube_PiecewiseFromRatio.fits',
                              'TemperatureCube_DendrogramObjects{0}_Piecewise.fits'),
                             ('','dendro')):

    for smooth in ("","_smooth",):#"_vsmooth"):
        cube = SpectralCube.read(hpath(ftemplate.format(smooth)))
        tproj = cube.mean(axis=1)
        tproj.wcs.wcs.ctype[0] = 'OFFSET'
        hdu = tproj.hdu
        hdu.header['CTYPE1'] = 'GLON'

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
