import pylab as pl
import glob
from astropy.io import fits
import aplpy
import paths
import os

if not os.path.isdir(paths.molpath('integ_imgs')):
    os.mkdir(paths.molpath('integ_imgs'))

fig1 = pl.figure(1)
for fn in glob.glob(paths.molpath("*integ.fits")):
    hdr = fits.getheader(fn)
    if hdr['CRVAL2'] == 0:
        fig1.clf()
        F = aplpy.FITSFigure(fn, convention='calabretta', figure=fig1)
        F.show_colorscale(cmap=pl.cm.hot)
        F.save(paths.molpath('integ_imgs/{0}.pdf').format(os.path.basename(fn[:-5])))
