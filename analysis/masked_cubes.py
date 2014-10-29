from spectral_cube import SpectralCube,BooleanArrayMask
from astropy import units as u
from paths import hpath
from astropy.io import fits
import time
from astropy import log

t0 = time.time()

cube303 = SpectralCube.read(hpath('APEX_H2CO_303_202_bl.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
cube321 = SpectralCube.read(hpath('APEX_H2CO_321_220_bl.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
mask = (fits.getdata(hpath('APEX_H2CO_303_202_mask.fits')).astype('bool') &
        cube303.mask.include(cube303._data, cube303.wcs) &
        cube321.mask.include(cube321._data, cube321.wcs))
bmask = BooleanArrayMask(mask, cube303.wcs)
cube303m = cube303.with_mask(bmask)
cube321m = cube321.with_mask(bmask)

cube303sm = SpectralCube.read(hpath('APEX_H2CO_303_202_smooth_bl.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
cube321sm = SpectralCube.read(hpath('APEX_H2CO_321_220_smooth_bl.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
masksm = (fits.getdata(hpath('APEX_H2CO_303_202_smooth_bl_mask.fits')).astype('bool') &
          cube303sm.mask.include(cube303sm._data, cube303sm.wcs) &
          cube321sm.mask.include(cube321sm._data, cube321sm.wcs))
bmasksm = BooleanArrayMask(masksm, cube303sm.wcs)
cube303msm = cube303sm.with_mask(bmasksm)
cube321msm = cube321sm.with_mask(bmasksm)


sncube = SpectralCube.read(hpath('APEX_H2CO_303_202_signal_to_noise_cube.fits'))
sncube._wcs = cube303._wcs
sncube.mask._wcs = cube303._wcs
sncubesm = SpectralCube.read(hpath('APEX_H2CO_303_202_smooth_signal_to_noise_cube.fits'))
sncubesm._wcs = cube303sm._wcs
sncubesm.mask._wcs = cube303sm._wcs

log.info("Masked cube creation took {0:0.1f} seconds".format(time.time()-t0))
