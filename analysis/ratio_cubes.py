import numpy as np
from spectral_cube import SpectralCube,BooleanArrayMask
from masked_cubes import (cube303m,cube321m,cube303,cube321,mask, bmask,
                          cube303msm, cube321msm, masksm, bmasksm,
                          sncube, sncubesm)
from noise import noise, noise_cube, sm_noise_cube

#noise = fits.getdata(hpath('APEX_H2CO_303_202_noise.fits'))
#noise = cube303[:50].std(axis=0).value
#noise = fits.getdata(mpath('APEX_H2CO_merge_high_sub_noise.fits'))
#nhits = nhits = fits.getdata(paths.mpath('APEX_H2CO_merge_high_nhits.fits'))
#noise[nhits<20] = np.nan

noise_flat = noise_cube[mask]
var_flat = noise_flat**2

ratio303321 = cube321m.flattened().value / cube303m.flattened().value
eratio303321 = ((ratio303321**2 * (var_flat/cube303m.flattened().value**2 +
                                   var_flat/cube321m.flattened().value**2))**0.5)
ratioOK = ratio303321 > eratio303321*3

data = np.zeros(cube303m.shape, dtype='float32')*np.nan
mask[mask] = ratioOK
data[mask] = ratio303321[ratioOK]
ratiocube_303321 = SpectralCube(data,
                                mask=BooleanArrayMask(np.isfinite(data),
                                                      wcs=cube303m.wcs),
                                wcs=cube303m.wcs)


noise_flat_sm = sm_noise_cube[masksm]
var_flat_sm = noise_flat_sm**2

ratio303321sm = cube321msm.flattened().value / cube303msm.flattened().value
eratio303321sm = ((ratio303321sm**2 * (var_flat_sm/cube303msm.flattened().value**2 +
                                       var_flat_sm/cube321msm.flattened().value**2))**0.5)
ratioOKsm = ratio303321sm > eratio303321sm*4

datasm = np.zeros(cube303msm.shape, dtype='float32')*np.nan
masksm[masksm] = ratioOKsm
datasm[masksm] = ratio303321sm[ratioOKsm]
ratiocubesm_303321 = SpectralCube(datasm,
                                  mask=BooleanArrayMask(np.isfinite(datasm),
                                                        wcs=cube303msm.wcs),
                                  wcs=cube303msm.wcs)
