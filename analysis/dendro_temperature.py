import time
import warnings
import numpy as np

import pylab as pl
from spectral_cube import SpectralCube, BooleanArrayMask
import FITS_tools
from numpy.lib.stride_tricks import as_strided
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy import log
from astropy import table
from astropy.utils.console import ProgressBar
from astrodendro import Dendrogram,ppv_catalog

from paths import hpath,mpath
from constrain_parameters import paraH2COmodel
from masked_cubes import cube303m,cube321m,cube303msm,cube321msm
from masked_cubes import mask as cube_signal_mask
from co_cubes import cube13co, cube18co, cube13cosm, cubeco18sm
from noise import noise, noise_cube, sm_noise_cube
from higal_gridded import column_regridded
from common_constants import logabundance,elogabundance

warnings.simplefilter('once')

if 'mf' not in locals():
    mf = paraH2COmodel()

# For debugging, to make it faster
# biggest_tree = dend[89]

def get_root(structure):
    """ Identify the root of the tree """
    if structure.parent is None:
        return structure
    else:
        return get_root(structure.parent)

def measure_dendrogram_properties(dend=None, cube303=cube303m,
                                  cube321=cube321m, cube13co=cube13co,
                                  cube18co=cube18co, noise_cube=noise_cube,
                                  suffix=""):

    assert (cube321.shape == cube303.shape == noise_cube.shape ==
            cube13co.shape == cube18co.shape)

    metadata = {}
    metadata['data_unit'] = u.K
    metadata['spatial_scale'] =  7.2 * u.arcsec
    metadata['beam_major'] =  30 * u.arcsec
    metadata['beam_minor'] =  30 * u.arcsec
    metadata['wavelength'] =  218.22219*u.GHz
    metadata['velocity_scale'] = u.km/u.s
    metadata['wcs'] = cube303.wcs

    keys = ['cmin1sig_chi2',
            'density_chi2',
            'dmin1sig_chi2',
            'cmax1sig_chi2',
            'tmax1sig_chi2',
            'column_chi2',
            'temperature_chi2',
            'tmin1sig_chi2',
            'eratio303321',
            'dmax1sig_chi2',
            'ratio303321',
            'logh2column',
            'elogh2column',
            'logabundance',
            'elogabundance',
           ]
    obs_keys = [
            'Stot303',
            'Stot321',
            'Smean303',
            'Smean321',
            'npix',
            'e303',
            'e321',
            'r303321',
            'er303321',
            '13cosum',
            'c18osum',
            '13comean',
            'c18omean',
            'parent',
            'root',
    ]
    columns = {k:[] for k in (keys+obs_keys)}

    log.debug("Initializing dendrogram temperature fitting loop")

    # Prepare an array to hold the fitted temperatures
    tcubedata = np.empty(cube303.shape, dtype='float32')
    tcubedata[:] = np.nan

    # FORCE wcs to match
    cube13co._wcs = cube18co._wcs = cube303.wcs
    cube13co.mask._wcs = cube18co.mask._wcs = cube303.wcs

    #objects = biggest_tree.descendants
    objects = dend
    catalog = ppv_catalog(objects, metadata)
    pb = ProgressBar(len(objects))
    for ii,structure in enumerate(objects):
        dend_obj_mask = BooleanArrayMask(structure.get_mask(), wcs=cube303.wcs)
        view = cube303.subcube_slices_from_mask(dend_obj_mask)
        submask = dend_obj_mask[view]
        assert submask.include().sum() == dend_obj_mask.include().sum()

        c303 = cube303[view].with_mask(submask)
        c321 = cube321[view].with_mask(submask)
        co13sum = cube13co.with_mask(dend_obj_mask).sum().value
        co18sum = cube18co.with_mask(dend_obj_mask).sum().value
        if hasattr(co13sum,'__len__'):
            raise TypeError(".sum() applied to an array has yielded a non scalar.")

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

        #error = (noise_cube[view][submask.include()]).sum() / submask.include().sum()**0.5
        var = ((noise_cube[dend_obj_mask.include()]**2).sum() / npix**2)
        error = var**0.5
        if np.isnan(error):
            raise ValueError("error is nan: this is impossible by definition.")
        er321303 = (r321303**2 * (var/Smean303**2 + var/Smean321**2))**0.5

        columns['Stot303'].append(Stot303)
        columns['Stot321'].append(Stot321)
        columns['Smean303'].append(Smean303)
        columns['Smean321'].append(Smean321)
        columns['npix'].append(npix)
        columns['e303'].append(error)
        columns['e321'].append(error)
        columns['r303321'].append(r321303)
        columns['er303321'].append(er321303)
        columns['13cosum'].append(co13sum)
        columns['c18osum'].append(co18sum)
        columns['13comean'].append(co13sum/npix)
        columns['c18omean'].append(co18sum/npix)
        columns['parent'].append(structure.parent.idx if structure.parent else -1)
        columns['root'].append(get_root(structure))

        mask2d = dend_obj_mask.include().max(axis=0)[view[1:]]
        logh2column = np.log10(np.nanmean(column_regridded.data[view[1:]][mask2d]) * 1e22)
        if np.isnan(logh2column):
            log.info("Source #{0} has NaNs".format(ii))
            logh2column = 24
        elogh2column = elogabundance

        if r321303 < 0 or np.isnan(r321303):
            for k in columns:
                if k not in obs_keys:
                    columns[k].append(np.nan)
        else:
            mf.set_constraints(ratio303321=r321303, eratio303321=er321303,
                               #ratio321322=ratio2, eratio321322=eratio2,
                               logh2column=logh2column, elogh2column=elogh2column,
                               logabundance=logabundance, elogabundance=elogabundance,
                               taline303=Smean303, etaline303=error,
                               taline321=Smean321, etaline321=error,
                               linewidth=5)
            row_data = mf.get_parconstraints()
            row_data['ratio303321'] = r321303
            row_data['eratio303321'] = er321303

            for k in row_data:
                columns[k].append(row_data[k])

            tcubedata[dend_obj_mask.include()] = row_data['temperature_chi2']

        if len(set(len(c) for k,c in columns.iteritems())) != 1:
            print("Columns are different lengths.  This is not allowed.")
            import ipdb; ipdb.set_trace()

        if ii % 100 == 0 or ii < 50:
            try:
                log.info("T: [{tmin1sig_chi2:7.2f},{temperature_chi2:7.2f},{tmax1sig_chi2:7.2f}]"
                         "  R={ratio303321:8.4f}+/-{eratio303321:8.4f}"
                         "  Smean303={Smean303:8.4f} +/- {e303:8.4f}"
                         "  Stot303={Stot303:8.2e}  npix={npix:6d}"
                         .format(Smean303=Smean303, Stot303=Stot303,
                                 npix=npix, e303=error, **row_data))

                pl.clf()
                mf.denstemplot()
                pl.draw()
                pl.show()
            except:
                pass
        else:
            pb.update(ii+1)

    for k in columns:
        if k not in catalog.keys():
            catalog.add_column(table.Column(name=k, data=columns[k]))

    for mid,lo,hi in (('temperature_chi2','tmin1sig_chi2','tmax1sig_chi2'),
                      ('density_chi2','dmin1sig_chi2','dmax1sig_chi2'),
                      ('column_chi2','cmin1sig_chi2','cmax1sig_chi2')):
        catalog.add_column(table.Column(name='elo_'+mid[0],
                                        data=catalog[mid]-catalog[lo]))
        catalog.add_column(table.Column(name='ehi_'+mid[0],
                                        data=catalog[hi]-catalog[mid]))

    catalog.write(hpath('PPV_H2CO_Temperature{0}.ipac'.format(suffix)), format='ascii.ipac')

    # Note that there are overlaps in the catalog, which means that ORDER MATTERS
    # in the above loop.  I haven't yet checked whether large scale overwrites
    # small or vice-versa; it may be that both views of the data are interesting.
    tcube = SpectralCube(data=tcubedata, wcs=cube303.wcs,
                         mask=cube303.mask, meta={'unit':'K'},
                         header=cube303.header,
                        )

    tcube.write(hpath('TemperatureCube_DendrogramObjects{0}.fits'.format(suffix)),
                overwrite=True)

def do_dendro_temperatures_both():
    do_dendro_temperatures_sharp()
    do_dendro_temperatures_smooth()

def do_dendro_temperatures_sharp():
    from dendrograms import dend
    measure_dendrogram_properties(dend=dend, cube303=cube303m, cube321=cube321m,
                                  suffix="")

def do_dendro_temperatures_smooth():
    from dendrograms import dendsm
    measure_dendrogram_properties(dend=dendsm, cube303=cube303msm,
                                  cube321=cube321msm, cube13co=cube13cosm,
                                  cube18co=cube18cosm,
                                  noise_cube=sm_noise_cube,
                                  suffix="_smooth")
