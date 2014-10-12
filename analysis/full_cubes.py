from make_apex_cubes import all_lines
from spectral_cube import SpectralCube,BooleanArrayMask
import pyspeckit
from astropy.io import fits
from numpy.lib.stride_tricks import as_strided
from paths import mpath
import numpy as np
import time
from astropy import log
from astropy import constants
from astropy import units as u

cube_merge_high = SpectralCube.read(mpath('APEX_H2CO_merge_high_plait_all.fits'))
noise = fits.getdata(mpath('APEX_H2CO_merge_high_plait_all_noise.fits'))
nhits = fits.getdata(mpath('APEX_H2CO_merge_high_nhits.fits'))
noise[nhits<20] = np.nan
noise_cube = as_strided(noise, shape=cube_merge_high.shape,
                        strides=(0,)+noise.strides)

f1 = all_lines['H2CO_303_202']
f2 = all_lines['H2CO_321_220']
h2co_cube_merge_high = cube_merge_high.spectral_slab(f1*(1-(150*u.km/u.s/constants.c)),
                                                     f2*(1+(100*u.km/u.s/constants.c)))
h2co_noise_cube = noise_cube.spectral_slab(f1*(1-(150*u.km/u.s/constants.c)),
                                           f2*(1+(100*u.km/u.s/constants.c)))

pcube_merge_high = pyspeckit.Cube(cube=h2co_cube_merge_high._data,
                                  errorcube=h2co_noise_cube,
                                  header=h2co_cube_merge_high.header,
                                  xarr=h2co_cube_merge_high.spectral_axis,
                                 )
pcube_merge_high.xarr.refX = 218.22219
pcube_merge_high.xarr.refX_units = 'GHz'
