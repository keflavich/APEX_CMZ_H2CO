import time
import numpy as np
from numpy.polynomial.polynomial import Polynomial

from astropy import log
from spectral_cube import SpectralCube, BooleanArrayMask
from paths import hpath, apath
from astropy.table import Table, Column
from ratio_cubes import ratio303321, eratio303321, noise_flat
from masked_cubes import cube303m,cube321m,cube303,cube321,mask,cube303msm,cube321msm
from astropy.utils.console import ProgressBar
from astrodendro import Dendrogram,ppv_catalog
from astropy.io import fits
from dendrograms import dend, dendsm

fit_table = Table.read(apath('piecewise_tvsratio_fit.ipac'), format='ascii.ipac')

# Define a piecewise interpolated function...
pwtem = lambda x: np.piecewise(x,
                               [(x>r['MinBound']) & (x<r['MaxBound'])
                                for r in fit_table],
                               [Polynomial([r['const'], r['xcoef'], r['x2coef']])
                                for r in fit_table] + [lambda y: np.nan]
                               )

tcubedata = np.empty(cube303m.shape)
tcubedata[~mask] = np.nan
tcubedata[mask] = pwtem(ratio303321)

tcube = SpectralCube(data=tcubedata, wcs=cube303m.wcs,
                     mask=cube303m.mask, meta={'unit':'K'},
                     header=cube303m.header,
                    )
assert tcube.header['CUNIT3'] == 'km s-1'

tcube.write(hpath('TemperatureCube_PiecewiseFromRatio.fits'), overwrite=True)

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


# Try the same thing but on the dendrogrammed data.  Basically, this is 
# dendro_temperature but *much* faster

for sm,cubeA,cubeB,objects in zip(("","_smooth"),
                                  (cube303m,cube303msm),
                                  (cube321m,cube321msm),
                                  (dend,dendsm),):

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
