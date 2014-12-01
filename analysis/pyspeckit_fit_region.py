import numpy as np
from astropy import units as u
from pyspeckit_fitting import (simplemodel, simplemodel2, simple_fitter,
                               simple_fitter2, simple_fitter3)
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

regs = pyregion.open(paths.rpath('spectral_apertures.reg'))
regdict = {r.attr[1]['text']:r for r in regs}

def get_subregion_pcube(cube303m, cube303, cube321, region):
    #scube = cube_merge_high.subcube_from_ds9region(pyregion.ShapeList([region]))
    scube303m = cube303m.subcube_from_ds9region(pyregion.ShapeList([region]))
    scube303 = cube303.subcube_from_ds9region(pyregion.ShapeList([region]))
    scube321 = cube321.subcube_from_ds9region(pyregion.ShapeList([region]))
    # TODO: get error map
    #pcube = pyspeckit.Cube(cube=scube)
    pcube303 = pyspeckit.Cube(cube=scube303)
    pcube303.xarr.refX = cube303.wcs.wcs.restfrq
    pcube303.xarr.refX_units = 'Hz'
    pcube321 = pyspeckit.Cube(cube=scube321)
    pcube321.xarr.refX = cube321.wcs.wcs.restfrq
    pcube321.xarr.refX_units = 'Hz'
    pcube = pyspeckit.CubeStack([pcube303,pcube321,])
    pcube.specfit.Registry.add_fitter('h2co_simple', simple_fitter3, 4,
                                      multisingle='multi')
    pcube.xarr.refX = cube303m.wcs.wcs.restfrq
    pcube.xarr.refX_units = 'Hz'
    return pcube, scube303m

def do_g047_box():
    # That's awesome.
    pc,c3 = get_subregion_pcube(cube303, cube303, cube321, regdict['G0.47-0.07box'])
    c3slab = c3.spectral_slab(50*u.km/u.s, 125*u.km/u.s)
    moments = c3slab.moment1(axis=0)
    peak = c3slab.max(axis=0)
    vguesses = moments
    vguesses[vguesses.value<50] = np.nan
    vguesses[vguesses.value>125] = np.nan
    mask = peak.value>0.22
    vguesses[~mask] = np.nan
    do_pyspeck_fits_1comp(pc, m1=vguesses, vrange=(50,125),
                          limits=[(0,20), (50,125),(1,40),(0,1)],
                          guesses=[1,None,5,0.5],
                          start_from_point=(22,6))
    return pc


def do_pyspeck_fits_1comp(pcube, cube303m=None, vguesses='moment',
                          guesses=[1,None,5,0.5,0.7,1], vrange=(-105,125),
                          m1=None,
                          limits=[(0,20), (-105,125),(1,40),(0,1),(0.3,1.1),(0,1e5)],
                          **kwargs):

    if m1 is None:
        m1 = cube303m.moment1(axis=0).to(u.km/u.s)
    if vguesses == 'moment':
        guesses_simple = np.array([guesses[0:1]+[x]+guesses[2:] for x in m1.value.flat]).T.reshape((len(guesses),)+m1.shape)
        bad = (guesses_simple[1,:,:] < vrange[0]) | (guesses_simple[1,:,:] > vrange[1])
        guesses_simple[1,bad] = 25
    else:
        guesses_simple = guesses

    # Need to show pcube which pixels to ignore: m1 has nans
    pcube.mapplot.plane = m1.value
    pcube.fiteach(fittype='h2co_simple', multifit=True,
                  guesses=guesses_simple,
                  limited=[(True,True)] * len(guesses_simple),
                  limits=limits,
                  multicore=8,
                  integral=False,
                  **kwargs
                 )

    pcube.mapplot(estimator=0)
    pcube_orig = pcube.parcube.copy()
    pcube2 = remove_bad_pars(pcube.parcube, pcube.errcube, len(guesses),
                             min_nsig=4)

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
