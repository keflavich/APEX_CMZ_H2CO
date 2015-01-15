import time
import numpy as np
from numpy.polynomial.polynomial import Polynomial

from astropy import log
from spectral_cube import SpectralCube, BooleanArrayMask
from paths import hpath, apath
from astropy.table import Table, Column
from ratio_cubes import (ratio303321, eratio303321, noise_flat, mask, ratioOK,
                         ratio303321sm, ratioOKsm, masksm)
from masked_cubes import (cube303m,cube321m,cube303msm,cube321msm,
                          cube303,cube321,cube303sm,cube321sm,
                          sncube, sncubesm)
from astropy.utils.console import ProgressBar
from astrodendro import Dendrogram,ppv_catalog
from astropy.io import fits
from dendrograms import dend, dendsm, dend321, dend321sm
from piecewise_rtotem import pwtem
from temperature_cubes import tcubesm_direct, tcube_direct


tcube_direct.write(hpath('TemperatureCube_PiecewiseFromRatio.fits'), overwrite=True)

tcubesm_direct.write(hpath('TemperatureCube_smooth_PiecewiseFromRatio.fits'), overwrite=True)


# Do the integrated version
for smooth in ('', '_smooth'):
    top = hpath('APEX_H2CO_321_220{0}_bl_mask_integ.fits'.format(smooth))
    bot = hpath('APEX_H2CO_303_202{0}_bl_mask_integ.fits'.format(smooth))
    ff = fits.open(top)
    rr = ff[0].data/fits.getdata(bot)
    tem = pwtem(rr.ravel())
    ff[0].data = tem.reshape(rr.shape)
    ff.writeto(hpath('IntegratedRatioPiecewiseTemperature{0}_303to321.fits'.format(smooth)),
               clobber=True)

# Moved from temperature_cube: make integrated temperature maps
for sm in ("","_smooth",'_321','_321smooth'):
    outpath = 'TemperatureCube_DendrogramObjects{0}_Piecewise.fits'.format(sm)
    tcube = SpectralCube.read(hpath(outpath))

    integ = tcube.mean(axis=0)
    integ.hdu.writeto(hpath(outpath).replace(".fits","_integ.fits"),
                      clobber=True)

hdu_template = integ.hdu

# integrated temperature maps weighted by 303 brightness
for sm in ("","_smooth",'_321','_321smooth'):
    outpath = 'TemperatureCube_DendrogramObjects{0}_Piecewise.fits'.format(sm)
    weight_cube = cube303msm if 'smooth' in sm else cube303m
    weights = weight_cube.filled_data[:].value
    tcube = fits.getdata(hpath(outpath))

    integ = np.nansum(tcube*weights, axis=0) / np.nansum(weights, axis=0)
    hdu_template.data = integ

    hdu_template.writeto(hpath(outpath).replace(".fits","_integ_weighted.fits"),
                         clobber=True)


# Try the same thing but on the dendrogrammed data.  Basically, this is 
# dendro_temperature but *much* faster

for sm,cubeA,cubeB,objects in zip(("","_smooth",'_321','_321smooth'),
                                  (cube303,cube303sm,cube303,cube303sm),
                                  (cube321,cube321sm,cube321,cube321sm),
                                  (dend,dendsm,dend321,dend321sm),):

    # reset data
    tcubedata = np.empty(cubeA.shape)
    tcubedata[:] = np.nan
    rcubedata = np.empty(cubeA.shape)
    rcubedata[:] = np.nan

    pb = ProgressBar(len(objects))
    for ii,structure in enumerate(objects):
        dend_obj_mask = BooleanArrayMask(structure.get_mask(), wcs=cubeA.wcs)
        view = cubeA.subcube_slices_from_mask(dend_obj_mask)
        submask = dend_obj_mask[view]
        assert submask.include().sum() == dend_obj_mask.include().sum()

        c303 = cubeA[view].with_mask(submask)
        c321 = cubeB[view].with_mask(submask)

        npix = submask.include().sum()
        Stot303 = c303.sum().value
        Stot321 = c321.sum().value
        if npix == 0:
            raise ValueError("npix=0. This is impossible.")
        Smean303 = Stot303/npix
        Smean321 = Stot321/npix
        try:
            r321303 = Stot321/Stot303
        except ZeroDivisionError:
            # py.FuckOff...
            r321303 = np.nan

        rcubedata[structure.get_mask()] = r321303
        tcubedata[structure.get_mask()] = pwtem(np.array([r321303]))

        pb.update(ii+1)

    # Note that there are overlaps in the catalog, which means that ORDER MATTERS
    # in the above loop.  I haven't yet checked whether large scale overwrites
    # small or vice-versa; it may be that both views of the data are interesting.
    tcube = SpectralCube(data=tcubedata, wcs=cubeA.wcs,
                         mask=cubeA.mask, meta={'unit':'K'},
                         header=cubeA.header,
                        )

    outpath = 'TemperatureCube_DendrogramObjects{0}_Piecewise.fits'.format(sm)
    tcube.write(hpath(outpath), overwrite=True)

    rcube = SpectralCube(data=rcubedata, wcs=cubeA.wcs,
                         mask=cubeA.mask, meta={'unit':'K'},
                         header=cubeA.header,
                        )

    outpath = 'RatioCube_DendrogramObjects{0}.fits'.format(sm)
    rcube.write(hpath(outpath), overwrite=True)

    max_temcube = tcube.max(axis=0)
    max_temcube.hdu.writeto(hpath('TemperatureCube_DendrogramObjects{0}_Piecewise_max.fits'.format(sm)), clobber=True)
    max_rcube = rcube.max(axis=0)
    max_rcube.hdu.writeto(hpath('RatioCube_DendrogramObjects{0}_Piecewise_max.fits'.format(sm)), clobber=True)

    mean_temcube = tcube.mean(axis=0)
    mean_temcube.hdu.writeto(hpath('TemperatureCube_DendrogramObjects{0}_Piecewise_mean.fits'.format(sm)), clobber=True)
    mean_rcube = rcube.mean(axis=0)
    mean_rcube.hdu.writeto(hpath('RatioCube_DendrogramObjects{0}_Piecewise_mean.fits'.format(sm)), clobber=True)

    hdu_template = mean_rcube.hdu
    tcube = tcube.filled_data[:].value
    weight_cube = cube303sm if 'smooth' in sm else cube303
    weights = weight_cube.filled_data[:].value
    weights[weights < 0] = 0

    mean_temcube = np.nansum(tcube*weights, axis=0) / np.nansum(weights, axis=0)
    hdu_template.data = mean_temcube
    hdu_template.writeto(hpath('TemperatureCube_DendrogramObjects{0}_Piecewise_weightedmean.fits'.format(sm)), clobber=True)
    mean_rcube = np.nansum(rcube*weights,axis=0) / np.nansum(weights, axis=0)
    hdu_template.data = mean_rcube
    hdu_template.writeto(hpath('RatioCube_DendrogramObjects{0}_Piecewise_weightedmean.fits'.format(sm)), clobber=True)
