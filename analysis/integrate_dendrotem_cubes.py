import numpy as np
from paths import hpath,mpath,fpath
from astropy import log
from spectral_cube import SpectralCube, BooleanArrayMask
from masked_cubes import (cube303m,cube321m,cube303msm,cube321msm,
                          cube303,cube321,cube303sm,cube321sm,
                          sncube, sncubesm)


for suffix in ("","_smooth"):

    outpath = 'TemperatureCube_DendrogramObjects{0}.fits'
    tcube = SpectralCube.read(hpath(outpath.format(suffix)))

    outpath_leaf = 'TemperatureCube_DendrogramObjects{0}_leaves.fits'
    tcubeleaf = SpectralCube.read(hpath(outpath_leaf.format(suffix)))

    integ = tcube.mean(axis=0)
    integ.hdu.writeto(hpath(outpath.format(suffix)).replace(".fits","_integ.fits"),
                      clobber=True)
    integleaf = tcubeleaf.mean(axis=0)
    integleaf.hdu.writeto(hpath(outpath_leaf.format(suffix)).replace(".fits","_integ.fits"),
                          clobber=True)

    hdu_template = integ.hdu

    log.info("Writing Weighted Integrated TemperatureCube")
    tcubed = tcube.filled_data[:].value
    weight_cube = cube303sm if 'smooth' in suffix else cube303
    weights = weight_cube.filled_data[:].value
    weights[weights < 0] = 0
    weights[np.isnan(tcubed)] = 0
    assert tcubed.shape == weights.shape
    mean_tem = np.nansum(tcubed*weights,axis=0) / np.nansum(weights, axis=0)
    hdu_template.data = mean_tem
    hdu_template.writeto(hpath(outpath.format(suffix)).replace(".fits","_integ_weighted.fits"),
                         clobber=True)

    log.info("Writing Weighted Integrated TemperatureCube (leaves only)")
    tcubedleaf = tcubeleaf.filled_data[:].value
    weights[np.isnan(tcubedleaf)] = 0
    mean_tem_leaf = np.nansum(tcubedleaf*weights,axis=0) / np.nansum(weights, axis=0)
    hdu_template.data = mean_tem_leaf
    hdu_template.writeto(hpath(outpath_leaf.format(suffix)).replace(".fits","_integ_weighted.fits"),
                         clobber=True)
