import numpy as np
from numpy.polynomial.polynomial import Polynomial

from astropy import log
from spectral_cube import SpectralCube, BooleanArrayMask
from paths import hpath, apath
from astropy.table import Table, Column
from ratio_cubes import ratio303321, eratio303321, noise_flat
from masked_cubes import cube303m,cube321m,cube303,cube321,mask

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

tcube.write(hpath('TemperatureCube_PiecewiseFromRatio.fits'), overwrite=True)

log.warn("This approach is far too noisy.  Need to do "
         "something like Voronoi tesselation or at least "
         "smoothing/downsampling to get adequate S/N.")
