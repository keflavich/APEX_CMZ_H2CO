from spectral_cube import SpectralCube,BooleanArrayMask
import pyspeckit
from astropy.io import fits
from numpy.lib.stride_tricks import as_strided
from paths import mpath
import numpy as np
import time
from astropy import log

cube_merge_high = SpectralCube.read(mpath('APEX_H2CO_merge_high_plait_all.fits'))
noise = fits.getdata(mpath('APEX_H2CO_merge_high_plait_all_noise.fits'))
nhits = fits.getdata(mpath('APEX_H2CO_merge_high_nhits.fits'))
noise[nhits<20] = np.nan
noise_cube = as_strided(noise, shape=cube_merge_high.shape,
                        strides=(0,)+noise.strides)

pcube_merge_high = pyspeckit.Cube(cube=cube_merge_high._data,
                                  errorcube=noise_cube,
                                  header=cube_merge_high.header,
                                  xarr=cube_merge_high.spectral_axis,
                                 )
pcube_merge_high.xarr.refX = 218.22219
pcube_merge_high.xarr.refX_units = 'GHz'
