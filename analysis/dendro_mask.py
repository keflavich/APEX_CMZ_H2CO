"""
Create an H2CO cloud catalog using dendrograms for further analysis using the temperature fitter
"""
import time

import numpy as np
from astropy import log
from astropy.io import fits

from paths import hpath,mpath
from astrodendro import Dendrogram
from noise import noise, noise_cube, sm_noise
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
              outfn="DendroMask_H2CO303202.hdf5"):

    # Use a little sigma-rejection to get a decently robust noise estimate
    noise_std = noise[noise==noise].std()
    noise_mean = noise[noise==noise].mean()
    err_estimate = noise[(noise > (noise_mean-noise_std)) &
                         (noise < (noise_mean+noise_std))].std()

    dend = Dendrogram.compute(cube, min_value=3*err_estimate,
                              min_delta=2*err_estimate, min_npix=50,
                              verbose=True)

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
    dend = make_dend(masked_cubes.cube303, noise)
    t2 = time.time()
    log.info("Dendrogramming took {0:0.1f} seconds".format(t2-t1))

    #t0 = time.time()
    #sncube_sm = make_sncube(smooth=True)
    t1 = time.time()
    #log.info("Smooth S/N cubemaking took {0:0.1f} seconds".format(t1-t0))
    dendsm = make_dend(masked_cubes.cube303sm, sm_noise,
                       outfn="DendroMask_H2CO303202_smooth.hdf5")
    t2 = time.time()
    log.info("Smooth Dendrogramming took {0:0.1f} seconds".format(t2-t1))
