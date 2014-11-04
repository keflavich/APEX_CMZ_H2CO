"""
Create an H2CO cloud catalog using dendrograms for further analysis using the temperature fitter
"""
import time

import numpy as np
from astropy import log
from astropy.io import fits

from paths import hpath,mpath
from astrodendro import Dendrogram
from noise import noise, noise_cube, sm_noise, cube303nm, cube303nmsm
import masked_cubes

def make_sncube(write=True, smooth=False):

    suff = 'smooth_' if smooth else ''
    fn = 'APEX_H2CO_303_202_{0}bl.fits'.format(suff)
    n = sm_noise if smooth else noise
    sncube = fits.getdata(hpath(fn)) / n
    ff = fits.open(hpath(fn))
    ff[0].data = sncube

    if write:
        outfn = 'APEX_H2CO_303_202_{0}signal_to_noise_cube.fits'.format(suff)
        ff.writeto(hpath(outfn), clobber=True)

    return sncube

def make_sn_dend(sncube, view=True, write=True,
              outfn="DendroMask_H2CO303202_signal_to_noise.hdf5"):

    dend = Dendrogram.compute(sncube, min_value=3, min_delta=2, min_npix=50,
                              verbose=True)

    if view:
        dend.viewer

    if write:
        dend.save_to(hpath(outfn))

    return dend

def make_dend(cube, noise, view=True, write=True,
              min_npix=100,
              outfn="DendroMask_H2CO303202.hdf5"):

    # Use a little sigma-rejection to get a decently robust noise estimate
    noise_std = noise[noise==noise].std()
    noise_mean = noise[noise==noise].mean()
    err_estimate = noise[(noise > (noise_mean-noise_std)) &
                         (noise < (noise_mean+noise_std))].mean()

    dend = Dendrogram.compute(cube.filled_data[:].value,
                              min_value=3*err_estimate,
                              min_delta=2*err_estimate,
                              min_npix=min_npix,
                              verbose=True, wcs=cube.wcs)

    if view:
        dend.viewer

    if write:
        dend.save_to(hpath(outfn))

    return dend

if __name__ == "__main__":
    #t0 = time.time()
    #sncube = make_sncube()
    t1 = time.time()
    #log.info("S/N cubemaking took {0:0.1f} seconds".format(t1-t0))
    dend = make_dend(cube303nm, noise)
    t2 = time.time()
    log.info("Dendrogramming {1} objects took {0:0.1f} seconds".format(t2-t1,
                                                                       len(dend)))

    #t0 = time.time()
    #sncube_sm = make_sncube(smooth=True)
    t1 = time.time()
    #log.info("Smooth S/N cubemaking took {0:0.1f} seconds".format(t1-t0))
    # I think the noise is higher than reported in the noise mask by some fraction;
    # there are "chunks" at the top of the map that are no good.
    dendsm = make_dend(cube303nmsm, sm_noise*1.5, min_npix=200,
                       outfn="DendroMask_H2CO303202_smooth.hdf5")
    t2 = time.time()
    log.info("Smooth Dendrogramming {1} objects"
             " took {0:0.1f} seconds".format(t2-t1, len(dendsm)))
