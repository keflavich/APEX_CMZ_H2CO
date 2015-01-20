"""
Piecewise interpolator from ratio 303_202 / 321_220 to kinetic temperature
(inferred from a set of dendrogram analyses)
"""
from astropy.table import Table, Column
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from paths import apath
from temperature_mapper import tm

fit_table = Table.read(apath('piecewise_tvsratio_fit.ipac'), format='ascii.ipac')

# old, wrong version # Define a piecewise interpolated function...
# old, wrong version pwtem = lambda x: np.piecewise(x,
# old, wrong version                                [(x>r['MinBound']) & (x<r['MaxBound'])
# old, wrong version                                 for r in fit_table],
# old, wrong version                                [Polynomial([r['const'], r['xcoef'], r['x2coef']])
# old, wrong version                                 for r in fit_table] + [lambda y: np.nan]

# Grabbed from dendrotem_plot on Jan 15, 2015
#pwtem = lambda x: np.polyval([190.12665966, 276.45454806,  11.09564855], x)
pwtem = lambda x: tm(x, density=10**4.5)
