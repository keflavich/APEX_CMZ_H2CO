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

ntemp,ndens,ncol = 5,30,30

temperatures = np.array([20,30,50,100,200])
densities = np.linspace(4,8,ndens)
columns = np.linspace(13, 18, ncol)
abundance = 1.2e-9 # Johnston / Ao
opr = 0.01 # assume primarily para
opr = 3
fortho = opr/(1+opr)

# execution time ~30s for a 5x10x10
TI,pars,bad_pars = compute_grid(Radex=pyradex.fjdu.Fjdu,
                                temperatures=temperatures,
                                densities=densities,
                                columns=columns,
                                run_kwargs={})

density_label = 'Density $n(\mathrm{H}_2)$ [log cm$^{-3}$]'
column_label = 'p-H$_2$CO [log cm$^{-2}$/(km s$^{-1}$ pc)]'

for ii in range(ntemp):
    fig = pl.figure(ii+1)
    fig.clf()
    ax = fig.gca()
    # Temperature, Density, Column
    im = ax.imshow(pars['fluxgrid_321'][ii,:,:]/pars['fluxgrid_303'][ii,:,:],
                   extent=[ # left, right, bottom, top
                           columns.min(), columns.max(),
                           densities.min(), densities.max(),
                          ],
                   cmap=pl.cm.gray_r)
    con = ax.contour(columns, densities,
                     pars['fluxgrid_321'][ii,:,:]/pars['fluxgrid_303'][ii,:,:],
                     levels=[0.5], colors=['r'])
    cb = fig.colorbar(mappable=im)
    ax.set_xlabel(column_label, labelpad=20)
    ax.set_ylabel(density_label)
    ax.set_title("Temperature $T={0}$ K".format(temperatures[ii]))
    cb.set_label("$3_{2,1}-2_{2,0}$ / $3_{0,3}-2_{0,2}$", labelpad=20)
    fig.savefig(paths.fpath("radex/opticallythick_T{0}K.pdf".format(temperatures[ii])))
