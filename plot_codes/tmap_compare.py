"""
Compare the various temperature maps to ensure that they are identical
"""
import numpy as np
from spectral_cube import SpectralCube, BooleanArrayMask
from temperature_cubes import tcubesm_direct, tcube_direct
from astropy.io import fits
import paths
from paths import hpath
import copy

""" Generated in make_ratiotem_cubesims.py """
im1 = fits.getdata(hpath('H2CO_321220_to_303202_bl_integ_temperature_dens1e4_masked.fits'))
im2 = fits.getdata(hpath('H2CO_321220_to_303202_bl_integ_temperature_dens1e4.fits'))
im3 = fits.getdata(hpath('H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e4_masked.fits'))

tcube = tcube_direct
tcube_mean = tcube.mean(axis=0)

wcube = SpectralCube.read(hpath('APEX_H2CO_303_202_bl.fits'))
mcube = SpectralCube.read(hpath('APEX_H2CO_303_202_bl_mask.fits'))

if tcube.shape != wcube.shape:
    raise

weighted = copy.copy(tcube)
weighted._data = wcube._data * tcube._data
mask = (wcube.mask & tcube.mask &
        BooleanArrayMask(weighted.filled_data[...] != 0, weighted.wcs) &
        BooleanArrayMask(wcube.filled_data[...] != 0, wcube.wcs) &
        BooleanArrayMask(mcube._data==1, mcube.wcs)
       )
weighted = weighted.with_mask(mask)
wcube = wcube.with_mask(mask)

tcube_wmean = tcube.mean(axis=0)

pl.figure(14, figsize=(12,20)).clf()
pl.subplot(5,1,1).imshow(im1)
pl.subplot(5,1,2).imshow(im2)
pl.subplot(5,1,3).imshow(im3)
pl.subplot(5,1,4).imshow(tcube_mean.value)
pl.subplot(5,1,5).imshow(tcube_wmean.value)
