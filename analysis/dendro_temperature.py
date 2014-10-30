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
from masked_cubes import (cube303m,cube321m,cube303msm,cube321msm,
                          cube303,cube321,cube303sm,cube321sm,
                          sncube, sncubesm)
from masked_cubes import mask as cube_signal_mask
from co_cubes import cube13co, cube18co, cube13cosm, cube18cosm
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

def measure_dendrogram_properties(dend=None, cube303=cube303,
                                  cube321=cube321, cube13co=cube13co,
                                  cube18co=cube18co, noise_cube=noise_cube,
                                  sncube=sncube,
                                  suffix="",
                                  last_index=None,
                                  write=True):

    assert (cube321.shape == cube303.shape == noise_cube.shape ==
            cube13co.shape == cube18co.shape == sncube.shape)
    assert sncube.wcs is cube303.wcs is sncube.mask._wcs

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
            's_ntotal',
            'index',
            'parent',
            'root',
            'lon',
            'lat',
            'vcen',
    ]
    columns = {k:[] for k in (keys+obs_keys)}

    log.debug("Initializing dendrogram temperature fitting loop")

    # Prepare an array to hold the fitted temperatures
    tcubedata = np.empty(cube303.shape, dtype='float32')
    tcubedata[:] = np.nan

    # FORCE wcs to match
    cube13co._wcs = cube18co._wcs = cube303.wcs
    cube13co.mask._wcs = cube18co.mask._wcs = cube303.wcs

    catalog = ppv_catalog(dend, metadata)
    pb = ProgressBar(len(catalog))
    for ii,row in enumerate(catalog):
        structure = dend[row['_idx']]
        assert structure.idx == row['_idx'] == ii
        dend_obj_mask = BooleanArrayMask(structure.get_mask(), wcs=cube303.wcs)
        dend_inds = structure.indices()

        view = cube303.subcube_slices_from_mask(dend_obj_mask)
        submask = dend_obj_mask[view]
        assert submask.include().sum() == dend_obj_mask.include().sum()

        sn = sncube[view].with_mask(submask)
        sntot = sn.sum().value
        #np.testing.assert_almost_equal(sntot, structure.values().sum(), decimal=0)

        c303 = cube303[view].with_mask(submask)
        c321 = cube321[view].with_mask(submask)
        co13sum = cube13co.with_mask(dend_obj_mask).sum().value
        co18sum = cube18co.with_mask(dend_obj_mask).sum().value
        if hasattr(co13sum,'__len__'):
            raise TypeError(".sum() applied to an array has yielded a non scalar.")

        npix = submask.include().sum()
        assert npix == structure.get_npix()
        Stot303 = c303.sum().value
        if np.isnan(Stot303):
            raise ValueError("NaN in cube.  This can't happen: the data from "
                             "which the dendrogram was derived can't have "
                             "NaN pixels.")

        Stot321 = c321.sum().value
        if npix == 0:
            raise ValueError("npix=0. This is impossible.")
        Smean303 = Stot303/npix
        if Stot303 <= 0:
            raise ValueError("The 303 flux is <=0.  This isn't possible because "
                             "the dendrogram was derived from the 303 data with a "
                             "non-zero threshold.")
        if np.isnan(Stot321):
            raise ValueError("NaN in 321 line")
        Smean321 = Stot321/npix

        #error = (noise_cube[view][submask.include()]).sum() / submask.include().sum()**0.5
        var = ((noise_cube[dend_obj_mask.include()]**2).sum() / npix**2)
        error = var**0.5
        if np.isnan(error):
            raise ValueError("error is nan: this is impossible by definition.")

        if Stot321 < 0:
            r321303 = error / Smean303
            er321303 = (r321303**2 * (var/Smean303**2 + 1))**0.5
        else:
            r321303 = Stot321 / Stot303
            er321303 = (r321303**2 * (var/Smean303**2 + var/Smean321**2))**0.5


        for c in columns:
            assert len(columns[c]) == ii

        columns['index'].append(row['_idx'])
        columns['s_ntotal'].append(sntot)
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
        s303 = cube303._data[dend_inds]
        x,y,z = cube303.world[dend_inds]
        lon = ((z.value-(360*(z.value>180)))*s303).sum()/s303.sum()
        lat = (y*s303).sum()/s303.sum()
        vel = (x*s303).sum()/s303.sum()
        columns['lon'].append(lon.value)
        columns['lat'].append(lat.value)
        columns['vcen'].append(vel.value)

        mask2d = dend_obj_mask.include().max(axis=0)[view[1:]]
        logh2column = np.log10(np.nanmean(column_regridded.data[view[1:]][mask2d]) * 1e22)
        if np.isnan(logh2column):
            log.info("Source #{0} has NaNs".format(ii))
            logh2column = 24
        elogh2column = elogabundance

        if r321303 < 0 or np.isnan(r321303):
            raise ValueError("Ratio <0: This can't happen any more because "
                             "if either num/denom is <0, an exception is "
                             "raised earlier")
            #for k in columns:
            #    if k not in obs_keys:
            #        columns[k].append(np.nan)
        else:
            # Replace negatives for fitting
            if Smean321 <= 0:
                Smean321 = error
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

        for c in columns:
            assert len(columns[c]) == ii+1

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

        if last_index is not None and ii >= last_index:
            break

    if last_index is not None:
        catalog = catalog[:last_index+1]

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

    if write:
        catalog.write(hpath('PPV_H2CO_Temperature{0}.ipac'.format(suffix)), format='ascii.ipac')

    # Note that there are overlaps in the catalog, which means that ORDER MATTERS
    # in the above loop.  I haven't yet checked whether large scale overwrites
    # small or vice-versa; it may be that both views of the data are interesting.
    tcube = SpectralCube(data=tcubedata, wcs=cube303.wcs,
                         mask=cube303.mask, meta={'unit':'K'},
                         header=cube303.header,
                        )

    if write:
        tcube.write(hpath('TemperatureCube_DendrogramObjects{0}.fits'.format(suffix)),
                    overwrite=True)

    return catalog, tcube

def do_dendro_temperatures_both():
    do_dendro_temperatures_sharp()
    do_dendro_temperatures_smooth()

def do_dendro_temperatures_sharp():
    from dendrograms import dend
    measure_dendrogram_properties(dend=dend, cube303=cube303, cube321=cube321,
                                  sncube=sncube, suffix="")

def do_dendro_temperatures_smooth():
    from dendrograms import dendsm
    assert sncubesm._wcs is cube303sm._wcs
    measure_dendrogram_properties(dend=dendsm, cube303=cube303sm,
                                  cube321=cube321sm, cube13co=cube13cosm,
                                  cube18co=cube18cosm,
                                  noise_cube=sm_noise_cube,
                                  sncube=sncubesm,
                                  suffix="_smooth")
