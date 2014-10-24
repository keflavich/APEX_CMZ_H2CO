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
noise[nhits<5] = np.nan
noise_cube = as_strided(noise, shape=cube_merge_high.shape,
                        strides=(0,)+noise.strides)
noise_spcube = SpectralCube(data=noise_cube, wcs=cube_merge_high.wcs)

cube_merge_high_sm = SpectralCube.read(mpath('APEX_H2CO_merge_high_plait_all_smooth.fits'))
noise_sm = fits.getdata(mpath('APEX_H2CO_merge_high_plait_all_smooth_noise.fits'))
noise_cube_sm = as_strided(noise_sm, shape=cube_merge_high_sm.shape,
                           strides=(0,)+noise_sm.strides)
noise_spcube_sm = SpectralCube(data=noise_cube_sm, wcs=cube_merge_high_sm.wcs)

# Create a cutout of the cube covering the H2CO lines
# it's barely worth it; cuts off 10% of pixels
f1 = all_lines['H2CO_303_202']*u.GHz
f2 = all_lines['H2CO_321_220']*u.GHz
h2co_cube_merge_high = cube_merge_high.spectral_slab(f1*(1-(150*u.km/u.s/constants.c)),
                                                     f2*(1+(100*u.km/u.s/constants.c)))
h2co_noise_cube = noise_spcube.spectral_slab(f1*(1-(150*u.km/u.s/constants.c)),
                                           f2*(1+(100*u.km/u.s/constants.c)))

h2co_cube_merge_high_sm = cube_merge_high_sm.spectral_slab(f1*(1-(150*u.km/u.s/constants.c)),
                                                           f2*(1+(100*u.km/u.s/constants.c)))
h2co_noise_cube_sm = noise_spcube_sm.spectral_slab(f1*(1-(150*u.km/u.s/constants.c)),
                                                 f2*(1+(100*u.km/u.s/constants.c)))


# Pyspeckit cube made from spectralcube
pcube_merge_high = pyspeckit.Cube(cube=h2co_cube_merge_high._data,
                                  errorcube=h2co_noise_cube._data,
                                  header=h2co_cube_merge_high.header,
                                  xarr=h2co_cube_merge_high.spectral_axis,
                                 )
pcube_merge_high.xarr.refX = 218.22219
pcube_merge_high.xarr.refX_units = 'GHz'

pcube_merge_high_sm = pyspeckit.Cube(cube=h2co_cube_merge_high_sm._data,
                                  errorcube=h2co_noise_cube_sm._data,
                                  header=h2co_cube_merge_high_sm.header,
                                  xarr=h2co_cube_merge_high_sm.spectral_axis,
                                 )
pcube_merge_high_sm.xarr.refX = 218.22219
pcube_merge_high_sm.xarr.refX_units = 'GHz'
