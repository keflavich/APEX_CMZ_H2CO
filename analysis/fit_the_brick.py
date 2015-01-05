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


brick_slice = np.s_[:,131:154,710:740]

brick_cube = cube_merge_high[brick_slice]
brick_pcube = pyspeckit.Cube(cube=brick_cube)
brick_pcube.xarr.refX = brick_cube.wcs.wcs.restfrq

brick_pcube.specfit.Registry.add_fitter('h2co_simple', simple_fitter2, 6,
                                        multisingle='multi')

b303m = cube303m[brick_slice]
m1 = b303m.moment1(axis=0).to(u.km/u.s)

# For radial multiscale fitting, see multiscale_fit.py

def do_pyspeck_fits_1comp():
    #guesses_simple = [1,25,5,0.5,0.7,1]
    guesses_simple = np.array([[1,x.value,5,0.5,0.7,1] for x in m1.flat]).T.reshape((6,)+m1.shape)
    bad = (guesses_simple[1,:,:] < -105) | (guesses_simple[1,:,:] > 125)
    guesses_simple[1,bad] = 25
    # Need to show brick_pcube which pixels to ignore: m1 has nans
    brick_pcube.mapplot.plane = m1.value
    brick_pcube.fiteach(fittype='h2co_simple', multifit=True,
                  guesses=guesses_simple,
                  limited=[(True,True)] * 6,
                  limits=[(0,20),(-105,125),(1,40),(0,1),(0.3,1.1),(0,1e5)],
                  multicore=8,
                 )

    brick_pcube.mapplot(estimator=0)

def do_pyspeck_fits_2comp():
    guesses_simple = [1,15,5,0.5,0.7,1] + [1,36,5,0.5,0.7,1]
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

def dendrify(view=True):
    from dendro_mask import make_dend
    return make_dend(cube=cube303nm[brick_slice], noise=noise, write=False,
                     view=view, min_nsig_delta=1, min_npix=50)

def dend_spec(xx=12, yy=17, min_nsig=3, min_delta=1, min_npix=3):
    b303 = cube303[brick_slice]
    n303 = noise[brick_slice[1:]]
    dend = Dendrogram.compute(b303[:,yy,xx].value, min_value=n303[yy,xx]*min_nsig,
                              min_delta=n303[yy,xx]*min_delta, min_npix=min_npix,
                              verbose=True,)
    dend_guesses = [b303.spectral_axis[leaf.indices()].mean() for leaf in dend]

    import pylab as pl
    pl.plot(b303[:,yy,xx].value, drawstyle='steps-mid', color='k', linewidth=0.5)
    for ii in range(len(dend)):
        pl.plot(dend[ii].indices()[0],
                b303[:,yy,xx].value[dend[ii].indices()], '.')

    return dend_guesses,dend

def dend_guess_npeaks(cube, noise_2d, min_nsig=3, min_delta=1, min_npix=3):
    valid_indices = np.where(np.isfinite(noise_2d))

    guesses = {}

    for yy,xx in ProgressBar(zip(*valid_indices)):
        dend = Dendrogram.compute(cube[:, yy, xx].value,
                                  min_value=min_nsig*noise_2d[yy,xx],
                                  min_delta=min_delta*noise_2d[yy,xx],
                                  min_npix=min_npix,)
        guesses[(yy,xx)] = [x
                            for leaf in dend
                            for x in 
                            (cube.spectral_axis[leaf.indices()].mean().value,
                             cube.spectral_axis[leaf.indices()].std().value)
                           ]
    return guesses

def guesses_to_cube(guesses, shape):
    maxlen = max((len(x) for x in guesses.values()))
    gcube = np.empty((maxlen,)+shape) + np.nan
    for (yy,xx),gg in guesses.iteritems():
        gcube[:len(gg),yy,xx] = gg
    return gcube
