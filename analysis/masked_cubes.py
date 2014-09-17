from spectral_cube import SpectralCube,BooleanArrayMask
from astropy import units as u
from paths import hpath
from astropy.io import fits
import time
from astropy import log

t0 = time.time()

cube303 = SpectralCube.read(hpath('APEX_H2CO_303_202.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
cube321 = SpectralCube.read(hpath('APEX_H2CO_321_220.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
mask = (fits.getdata(hpath('APEX_H2CO_303_202_mask.fits')).astype('bool') &
        cube303.mask.include(cube303._data, cube303.wcs) &
        cube321.mask.include(cube321._data, cube321.wcs))
bmask = BooleanArrayMask(mask, cube303.wcs)
cube303m = cube303.with_mask(bmask)
cube321m = cube321.with_mask(bmask)

log.debug("Masked cube creation took {0:0.1f} seconds".format(time.time()-t0))
