"""
Compare temperature maps made with order-of-magnitude different column,
density

(obsolete- not using these any more, I guess?)
(also, figures were never saved)
"""

#from astropy.io import fits
#from paths import hpath
#import pylab as pl
#import numpy as np
#
#im1 = fits.getdata(hpath('H2CO_321220_to_303202_bl_integ_temperature_dens1e5_col5e22.fits'))
#im2 = fits.getdata(hpath('H2CO_321220_to_303202_bl_integ_temperature_dens1e4_col5e22.fits'))
#im3 = fits.getdata(hpath('H2CO_321220_to_303202_bl_integ_temperature_dens1e5_col3e23.fits'))
#
#ok = np.isfinite(im1) & np.isfinite(im2) & np.isfinite(im3)
#
#fig1 = pl.figure(1)
#fig1.clf()
#pl.plot([0,300],[0,300], 'k--')
#pl.plot(im1[ok], im2[ok], ',')
#pl.plot(im1[ok], im3[ok], ',')
#
#fig2 = pl.figure(2)
#fig2.clf()
#pl.plot(im1[ok], (im1[ok]-im2[ok])/im1[ok], ',', label='Lower density')
#pl.plot(im1[ok], (im1[ok]-im3[ok])/im1[ok], ',', label='Higher column')
#pl.legend(loc='best')
