"""
Piecewise interpolator from ratio 303_202 / 321_220 to kinetic temperature
(inferred from a set of dendrogram analyses)
"""
from astropy.table import Table, Column
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from paths import apath

fit_table = Table.read(apath('piecewise_tvsratio_fit.ipac'), format='ascii.ipac')

# Define a piecewise interpolated function...
pwtem = lambda x: np.piecewise(x,
                               [(x>r['MinBound']) & (x<r['MaxBound'])
                                for r in fit_table],
                               [Polynomial([r['const'], r['xcoef'], r['x2coef']])
                                for r in fit_table] + [lambda y: np.nan]
                               )
