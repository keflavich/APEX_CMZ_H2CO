from __future__ import print_function
from astropy import log
import numpy as np
import pyradex
import copy
import os
import pylab as pl
import paths
from astropy import table
from paths import analysispath
import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.utils.console import ProgressBar
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))
from pyradex_h2comm_grid import compute_grid

ntemp,ndens,ncol = 5,10,10

temperatures = np.array([20,30,50,100,200])
densities = np.linspace(2.5,8,ndens)
columns = np.linspace(13, 18, ncol)
abundance = 1.2e-9 # Johnston / Ao
opr = 0.01 # assume primarily para
opr = 3
fortho = opr/(1+opr)

TI,pars,bad_pars = compute_grid(Radex=pyradex.fjdu.Fjdu,
                                temperatures=temperatures,
                                densities=densities,
                                columns=columns)

