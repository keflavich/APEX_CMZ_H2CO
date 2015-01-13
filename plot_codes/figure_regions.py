import numpy as np
import pyregion
from astropy import units as u
from paths import h2copath,figurepath,hpath,rpath,fpath
import copy
import os
import aplpy
import pylab as pl
import matplotlib
from astropy.io import fits
from masked_cubes import cube303, cube303m, cube321m, cube303msm, cube321msm, sncube, sncubesm
from astropy import log
import wcsaxes
from aplpy.regions import ds9
from astropy.wcs import WCS
from agpy import asinh_norm

pl.ioff()

cm = matplotlib.cm.RdYlBu_r
cm.set_bad('#888888')

if not os.path.exists(hpath("APEX_H2CO_303_202_integrated.fits")):
    im = cube303.sum(axis=0).hdu
    im.writeto(hpath("APEX_H2CO_303_202_integrated.fits"))
else:
    im = fits.open(hpath("APEX_H2CO_303_202_integrated.fits"))[0]

wcs = WCS(im.header)

fig1 = pl.figure(1)
fig1.clf()
ax = wcsaxes.WCSAxesSubplot(fig1, 1,1,1, wcs=wcs)
fig1.add_axes(ax)

cm = pl.cm.gray_r
im.data -= np.nanmin(im.data) - 0.00001
vmin = np.percentile(im.data[np.isfinite(im.data)], 00.05)
vmax = np.percentile(im.data[np.isfinite(im.data)], 99.95)
ax.imshow(im.data, cmap=cm, vmin=vmin, vmax=vmax,
          norm=matplotlib.colors.LogNorm())

region_file = rpath('spectral_apertures.reg')
regions = pyregion.open(region_file)

for sh in regions:
    sh.attr[1]['color'] = 'darkgreen'

boxes = pyregion.ShapeList([sh for sh in regions if sh.name == 'box'])
circles = pyregion.ShapeList([sh for sh in regions if sh.name == 'circle'])

PC, TC = ds9(boxes, im.header)
PC.add_to_axes(ax)
#TC.add_to_axes(ax)
pl.savefig(fpath("regions/boxes_on_h2co.pdf"), bbox_inches='tight')

PC.set_visible(False)
TC.set_visible(False)
PC, TC = ds9(circles, im.header)
PC.add_to_axes(ax)
#TC.add_to_axes(ax)
pl.savefig(fpath("regions/circles_on_h2co.pdf"), bbox_inches='tight')

regions = pyregion.open(rpath('target_fields_8x8_gal.reg'))
for sh in regions:
    sh.attr[1]['color'] = 'darkgreen'
PC.set_visible(False)
TC.set_visible(False)
PC, TC = ds9(regions, im.header)
PC.add_to_axes(ax)
#TC.add_to_axes(ax)
pl.savefig(fpath("regions/square_fields_on_h2co.pdf"), bbox_inches='tight')
