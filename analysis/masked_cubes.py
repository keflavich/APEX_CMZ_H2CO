import numpy as np
from spectral_cube import SpectralCube,BooleanArrayMask
from astropy import units as u
from paths import hpath
from astropy.io import fits
import time
from astropy import log

t0 = time.time()

hc3n_regions = [{'v':(-101,55),
                 'x':(500,533),
                 'y':(108,133),},
                {'v':(-133,-70),
                 'x':(787,884),
                 'y':(87,120),}]

def mask_out_region(mask_array, cube, regions=hc3n_regions):

    for region in regions:
        z = [cube.closest_spectral_channel(v*u.km/u.s)
             for v in region['v']]
        view = [slice(*z),
                slice(*region['y']),
                slice(*region['x'])
               ]

        mask_array[view] = False

    return mask_array

cube303 = SpectralCube.read(hpath('APEX_H2CO_303_202_bl.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
cube321 = SpectralCube.read(hpath('APEX_H2CO_321_220_bl.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
maskarr = mask_out_region(fits.getdata(hpath('APEX_H2CO_303_202_bl_mask.fits')).astype('bool'), cube303)
mask = (maskarr &
        cube303.mask.include(cube303, cube303.wcs) &
        cube321.mask.include(cube321, cube321.wcs))
bmask = BooleanArrayMask(mask, cube303.wcs)
cube303m = cube303.with_mask(bmask)
cube321m = cube321.with_mask(bmask)

cube303sm = SpectralCube.read(hpath('APEX_H2CO_303_202_smooth_bl.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
cube321sm = SpectralCube.read(hpath('APEX_H2CO_321_220_smooth_bl.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
smmaskarr = mask_out_region(fits.getdata(hpath('APEX_H2CO_303_202_smooth_bl_mask.fits')).astype('bool'), cube303sm)
masksm = (smmaskarr &
          cube303sm.mask.include(cube303sm, cube303sm.wcs) &
          cube321sm.mask.include(cube321sm, cube321sm.wcs))
bmasksm = BooleanArrayMask(masksm, cube303sm.wcs)
cube303msm = cube303sm.with_mask(bmasksm)
cube321msm = cube321sm.with_mask(bmasksm)

# resample smoothed mask onto original grid
masksm_rs = np.zeros_like(mask, dtype='bool')
masksm_rs[::2,:,:] = masksm
masksm_rs[1::2,:,:] = masksm
bmasksm_rs = BooleanArrayMask(masksm_rs, cube303.wcs)

sncube = SpectralCube.read(hpath('APEX_H2CO_303_202_signal_to_noise_cube.fits'))
sncube._wcs = cube303._wcs
sncube.mask._wcs = cube303._wcs
sncubesm = SpectralCube.read(hpath('APEX_H2CO_303_202_smooth_signal_to_noise_cube.fits'))
sncubesm._wcs = cube303sm._wcs
sncubesm.mask._wcs = cube303sm._wcs

log.info("Masked cube creation took {0:0.1f} seconds".format(time.time()-t0))
