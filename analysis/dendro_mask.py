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

def make_dend(sncube, view=True, write=True,
              outfn="DendroMask_H2CO303202_signal_to_noise.hdf5"):

    dend = Dendrogram.compute(sncube, min_value=3, min_delta=2, min_npix=50,
                              verbose=True)

    if view:
        dend.viewer

    if write:
        dend.save_to(hpath(outfn))

    return dend

if __name__ == "__main__":
    t0 = time.time()
    sncube = make_sncube()
    t1 = time.time()
    log.info("S/N cubemaking took {0:0.1f} seconds".format(t1-t0))
    dend = make_dend(sncube)
    t2 = time.time()
    log.info("Dendrogramming took {0:0.1f} seconds".format(t2-t1))

    t0 = time.time()
    sncube = make_sncube(smooth=True)
    t1 = time.time()
    log.info("Smooth S/N cubemaking took {0:0.1f} seconds".format(t1-t0))
    dend = make_dend(sncube,
                     outfn="DendroMask_H2CO303202_smooth_signal_to_noise.hdf5")
    t2 = time.time()
    log.info("Smooth Dendrogramming took {0:0.1f} seconds".format(t2-t1))
