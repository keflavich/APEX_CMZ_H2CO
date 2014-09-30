"""
Create an H2CO cloud catalog using dendrograms for further analysis using the temperature fitter
"""
import time

import numpy as np
from astropy import log
from astropy.io import fits

from paths import hpath,mpath
from astrodendro import Dendrogram
from noise import noise, noise_cube

def make_sncube(write=True):

    sncube = fits.getdata(hpath('APEX_H2CO_303_202_bl.fits')) / noise
    ff = fits.open(hpath('APEX_H2CO_303_202_bl.fits'))
    ff[0].data = sncube

    if write:
        ff.writeto(hpath('APEX_H2CO_303_202_signal_to_noise_cube.fits'),
                   clobber=True)

    return sncube

def make_dend(sncube, view=True, write=True):

    dend = Dendrogram.compute(sncube, min_value=3, min_delta=2, min_npix=50,
                              verbose=True)

    if view:
        dend.viewer

    if write:
        dend.save_to(hpath("DendroMask_H2CO303202_signal_to_noise.hdf5"))

    return dend

if __name__ == "__main__":
    t0 = time.time()
    sncube = make_sncube()
    t1 = time.time()
    log.debug("S/N cubemaking took {0:0.1f} seconds".format(t1-t0))
    dend = make_dend(sncube)
    t2 = time.time()
    log.debug("Dendrogramming took {0:0.1f} seconds".format(t2-t1))
