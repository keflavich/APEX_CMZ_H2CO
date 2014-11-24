import numpy as np
from astropy import units as u
from pyspeckit_fitting import (simplemodel, simplemodel2, simple_fitter,
                               simple_fitter2)
from full_cubes import cube_merge_high
from masked_cubes import (cube303, cube303sm, cube303m, cube321m, cube303msm,
                          cube321msm, cube321, cube321sm)
from noise import (noise, noise_cube, sm_noise, cube303nm, cube303nmsm,
                   cube321nm, cube321nmsm)
import pyspeckit
from astrodendro import Dendrogram
from astropy.utils.console import ProgressBar
import paths
import pylab as pl
import pyregion

def get_subregion_pcube(cube_merge_high, cube303m, region):
    scube = cube_merge_high.subcube_from_ds9region(pyregion.ShapeList([region]))
    scube303m = cube303m.subcube_from_ds9region(pyregion.ShapeList([region]))
    pcube = pyspeckit.Cube(cube=scube)
    pcube.specfit.Registry.add_fitter('h2co_simple', simple_fitter2, 6,
                                      multisingle='multi')
    pcube.xarr.refX = cube303m.wcs.wcs.restfrq
    return pcube, scube303m

def do_pyspeck_fits_1comp(pcube, cube303m, vguesses='moment',
                          guesses=[1,None,5,0.5,0.7,1], vrange=(-105,125)):

    m1 = cube303m.moment1(axis=0).to(u.km/u.s)
    if vguesses == 'moment':
        guesses_simple = np.array([guesses[0:1]+[x]+guesses[2:] for x in m1.value.flat]).T.reshape((6,)+m1.shape)
        bad = (guesses_simple[1,:,:] < vrange[0]) | (guesses_simple[1,:,:] > vrange[1])
        guesses_simple[1,bad] = 25
    else:
        guesses_simple = guesses

    # Need to show pcube which pixels to ignore: m1 has nans
    pcube.mapplot.plane = m1.value
    pcube.fiteach(fittype='h2co_simple', multifit=True,
                  guesses=guesses_simple,
                  limited=[(True,True)] * 6,
                  limits=[(0,20), vrange,(1,40),(0,1),(0.3,1.1),(0,1e5)],
                  multicore=8,
                 )

    pcube.mapplot(estimator=0)

def do_pyspeck_fits_2comp(pcube, cube303m, vrange=(-105,125), guesses_simple =
                          [1,15,5,0.5,0.7,1] + [1,36,5,0.5,0.7,1]):

    m1 = cube303m.moment1(axis=0).to(u.km/u.s)
    # Need to show brick_pcube which pixels to ignore: m1 has nans
    brick_pcube.mapplot.plane = m1.value
    brick_pcube.fiteach(fittype='h2co_simple', multifit=True,
                        guesses=guesses_simple, limited=[(True,True)] * 6,
                        limits=[(0,20),(-105,125),(1,40),(0,1),(0.3,1.1),(0,1e5)],
                        multicore=8,
                 )

    brick_pcube.mapplot(estimator=0)
    pcube_orig = brick_pcube.parcube.copy()
    pcube2 = remove_bad_pars(brick_pcube.parcube, brick_pcube.errcube, 6,
                             min_nsig=4)

def remove_bad_pars(parcube, errcube, npars, min_nsig=3):
    """
    Excise bad fits from a parameter fit cube (from pyspeckit)
    """
    if parcube.shape[0] % npars != 0:
        raise ValueError("Invalid # of parameters {0},"
                         "it should divide evenly "
                         "into {1}".format(npars, parcube.shape[0]))
    assert errcube.shape == parcube.shape

    ncomp = parcube.shape[0]/npars
    for jj in range(ncomp):
        shift = jj*npars
        bad = ((parcube[0+shift,:,:] < min_nsig*errcube[0+shift,:,:]) |
               (errcube[3+shift,:,:] == 0))
        parcube[jj*npars:(jj+1)*npars,bad] = 0
    return parcube
