from masked_cubes import cube303m,cube321m
from astropy.io import fits
from numpy.lib.stride_tricks import as_strided

noise = fits.getdata(mpath('APEX_H2CO_merge_high_sub_noise.fits'))
nhits = nhits = fits.getdata(paths.mpath('APEX_H2CO_merge_high_nhits.fits'))
noise[nhits<20] = np.nan
noise_cube = as_strided(noise, shape=mask.shape, strides=(0,)+noise.shape)
