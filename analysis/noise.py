from masked_cubes import cube303m,cube321m
from astropy.io import fits
from numpy.lib.stride_tricks import as_strided
from paths import mpath
import numpy as np
import time
from astropy import log

t0 = time.time()

noise = fits.getdata(mpath('APEX_H2CO_merge_high_plait_all_noise.fits'))
nhits = fits.getdata(mpath('APEX_H2CO_merge_high_nhits.fits'))
noise[nhits<20] = np.nan
noise_cube = as_strided(noise, shape=cube303m.shape, strides=(0,)+noise.strides)

sm_noise = fits.getdata(mpath('APEX_H2CO_merge_high_plait_all_smooth_noise.fits'))
sm_noise_cube = as_strided(sm_noise, shape=cube303m.shape, strides=(0,)+sm_noise.strides)

log.debug("Noise creation took {0:0.1f} seconds".format(time.time()-t0))
