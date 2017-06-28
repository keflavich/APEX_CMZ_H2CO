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
from astropy import constants
from astropy import log
from astropy import table
from astropy.utils.console import ProgressBar
from astrodendro import Dendrogram,ppv_catalog

from paths import hpath,mpath,fpath,tpath
from constrain_parameters import paraH2COmodel
from masked_cubes import (cube303m,cube321m,cube303msm,cube321msm,
                          cube303,cube321,cube303sm,cube321sm,
                          sncube, sncubesm)
from masked_cubes import mask as cube_signal_mask
from co_cubes import cube13co, cube18co, cube13cosm, cube18cosm
from noise import noise, noise_cube, sm_noise_cube
from higal_gridded import column_regridded, dusttem_regridded
from common_constants import logabundance,elogabundance
import heating
from gaussian_correction import gaussian_correction

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
                                  plot_some=True,
                                  line='303',
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

    keys = [
            'density_chi2',
            'expected_density',
            'dmin1sig_chi2',
            'dmax1sig_chi2',
            'column_chi2',
            'expected_column',
            'cmin1sig_chi2',
            'cmax1sig_chi2',
            'temperature_chi2',
            'expected_temperature',
            'tmin1sig_chi2',
            'tmax1sig_chi2',
            'eratio321303',
            'ratio321303',
            'logh2column',
            'elogh2column',
            'logabundance',
            'elogabundance',
           ]
    obs_keys = [
            'Stot303',
            'Smin303',
            'Smax303',
            'Stot321',
            'Smean303',
            'Smean321',
            'npix',
            'e303',
            'e321',
            'r321303',
            'er321303',
            '13cosum',
            'c18osum',
            '13comean',
            'c18omean',
            's_ntotal',
            'index',
            'is_leaf',
            'parent',
            'root',
            'lon',
            'lat',
            'vcen',
            'higaldusttem',
            'reff',
            'dustmass',
            'dustmindens',
            'bad',
            #'tkin_turb',
    ]
    columns = {k:[] for k in (keys+obs_keys)}

    log.debug("Initializing dendrogram temperature fitting loop")

    # FORCE wcs to match
    # (technically should reproject here)
    cube13co._wcs = cube18co._wcs = cube303.wcs
    cube13co.mask._wcs = cube18co.mask._wcs = cube303.wcs

    if line == '303':
        maincube = cube303
    elif line == '321':
        maincube = cube321
    else:
        raise ValueError("Unrecognized line: {0}".format(line))

    # Prepare an array to hold the fitted temperatures
    tcubedata = np.empty(maincube.shape, dtype='float32')
    tcubedata[:] = np.nan
    tcubeleafdata = np.empty(maincube.shape, dtype='float32')
    tcubeleafdata[:] = np.nan


    nbad = 0

    catalog = ppv_catalog(dend, metadata)
    pb = ProgressBar(len(catalog))
    for ii,row in enumerate(catalog):
        structure = dend[row['_idx']]
        assert structure.idx == row['_idx'] == ii
        dend_obj_mask = BooleanArrayMask(structure.get_mask(), wcs=cube303.wcs)
        dend_inds = structure.indices()

        view = (slice(dend_inds[0].min(), dend_inds[0].max()+1),
                slice(dend_inds[1].min(), dend_inds[1].max()+1),
                slice(dend_inds[2].min(), dend_inds[2].max()+1),)
        #view2 = cube303.subcube_slices_from_mask(dend_obj_mask)
        submask = dend_obj_mask[view]
        #assert np.count_nonzero(submask.include()) == np.count_nonzero(dend_obj_mask.include())

        sn = sncube[view].with_mask(submask)
        sntot = sn.sum().value
        #np.testing.assert_almost_equal(sntot, structure.values().sum(), decimal=0)

        c303 = cube303[view].with_mask(submask)
        c321 = cube321[view].with_mask(submask)
        co13sum = cube13co[view].with_mask(submask).sum().value
        co18sum = cube18co[view].with_mask(submask).sum().value
        if hasattr(co13sum,'__len__'):
            raise TypeError(".sum() applied to an array has yielded a non scalar.")

        npix = submask.include().sum()
        assert npix == structure.get_npix()
        Stot303 = c303.sum().value
        if np.isnan(Stot303):
            raise ValueError("NaN in cube.  This can't happen: the data from "
                             "which the dendrogram was derived can't have "
                             "NaN pixels.")
        Smax303 = c303.max().value
        Smin303 = c303.min().value

        Stot321 = c321.sum().value
        if npix == 0:
            raise ValueError("npix=0. This is impossible.")
        Smean303 = Stot303/npix
        if Stot303 <= 0 and line=='303':
            raise ValueError("The 303 flux is <=0.  This isn't possible because "
                             "the dendrogram was derived from the 303 data with a "
                             "non-zero threshold.")
        elif Stot303 <= 0 and line=='321':
            Stot303 = 0
            Smean303 = 0
        elif Stot321 <= 0 and line=='321':
            raise ValueError("The 321 flux is <=0.  This isn't possible because "
                             "the dendrogram was derived from the 321 data with a "
                             "non-zero threshold.")
        if np.isnan(Stot321):
            raise ValueError("NaN in 321 line")
        Smean321 = Stot321/npix

        #error = (noise_cube[view][submask.include()]).sum() / submask.include().sum()**0.5
        var = ((noise_cube[dend_obj_mask.include()]**2).sum() / npix**2)
        error = var**0.5
        if np.isnan(error):
            raise ValueError("error is nan: this is impossible by definition.")

        if line == '321' and Stot303 == 0:
            r321303 = np.nan
            er321303 = np.nan
        elif Stot321 < 0:
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
        columns['Smax303'].append(Smax303)
        columns['Smin303'].append(Smin303)
        columns['Stot321'].append(Stot321)
        columns['Smean303'].append(Smean303)
        columns['Smean321'].append(Smean321)
        columns['npix'].append(npix)
        columns['e303'].append(error)
        columns['e321'].append(error)
        columns['r321303'].append(r321303)
        columns['er321303'].append(er321303)
        columns['13cosum'].append(co13sum)
        columns['c18osum'].append(co18sum)
        columns['13comean'].append(co13sum/npix)
        columns['c18omean'].append(co18sum/npix)
        columns['is_leaf'].append(structure.is_leaf)
        columns['parent'].append(structure.parent.idx if structure.parent else -1)
        columns['root'].append(get_root(structure).idx)
        s_main = maincube._data[dend_inds]
        x,y,z = maincube.world[dend_inds]
        lon = ((z.value-(360*(z.value>180)))*s_main).sum()/s_main.sum()
        lat = (y*s_main).sum()/s_main.sum()
        vel = (x*s_main).sum()/s_main.sum()
        columns['lon'].append(lon)
        columns['lat'].append(lat.value)
        columns['vcen'].append(vel.value)

        mask2d = dend_obj_mask.include().max(axis=0)[view[1:]]
        logh2column = np.log10(np.nanmean(column_regridded.data[view[1:]][mask2d]) * 1e22)
        if np.isnan(logh2column):
            log.info("Source #{0} has NaNs".format(ii))
            logh2column = 24
        elogh2column = elogabundance
        columns['higaldusttem'].append(np.nanmean(dusttem_regridded.data[view[1:]][mask2d]))

        r_arcsec = row['radius']*u.arcsec
        reff = (r_arcsec*(8.5*u.kpc)).to(u.pc, u.dimensionless_angles())
        mass = ((10**logh2column*u.cm**-2)*np.pi*reff**2*2.8*constants.m_p).to(u.M_sun)
        density = (mass/(4/3.*np.pi*reff**3)/constants.m_p/2.8).to(u.cm**-3)

        columns['reff'].append(reff.value)
        columns['dustmass'].append(mass.value)
        columns['dustmindens'].append(density.value)
        mindens = np.log10(density.value)
        if mindens < 3:
            mindens = 3

        if (r321303 < 0 or np.isnan(r321303)) and line != '321':
            raise ValueError("Ratio <0: This can't happen any more because "
                             "if either num/denom is <0, an exception is "
                             "raised earlier")
            #for k in columns:
            #    if k not in obs_keys:
            #        columns[k].append(np.nan)
        elif (r321303 < 0 or np.isnan(r321303)) and line == '321':
            for k in keys:
                columns[k].append(np.nan)
        else:
            # Replace negatives for fitting
            if Smean321 <= 0:
                Smean321 = error
            mf.set_constraints(ratio321303=r321303, eratio321303=er321303,
                               #ratio321322=ratio2, eratio321322=eratio2,
                               logh2column=logh2column, elogh2column=elogh2column,
                               logabundance=logabundance, elogabundance=elogabundance,
                               taline303=Smean303, etaline303=error,
                               taline321=Smean321, etaline321=error,
                               mindens=mindens,
                               linewidth=10)
            row_data = mf.get_parconstraints()
            row_data['ratio321303'] = r321303
            row_data['eratio321303'] = er321303

            for k in row_data:
                columns[k].append(row_data[k])

            # Exclude bad velocities from cubes
            if row['v_cen'] < -80e3 or row['v_cen'] > 180e3:
                # Skip: there is no real structure down here
                nbad += 1
                is_bad = True
            else:
                is_bad = False
                tcubedata[dend_obj_mask.include()] = row_data['expected_temperature']
                if structure.is_leaf:
                    tcubeleafdata[dend_obj_mask.include()] = row_data['expected_temperature']

            columns['bad'].append(is_bad)

            width = row['v_rms']*u.km/u.s
            lengthscale = reff

            #REMOVED in favor of despotic version done in dendrograms.py
            # we use the analytic version here; the despotic version is
            # computed elsewhere (with appropriate gcor factors)
            #columns['tkin_turb'].append(heating.tkin_all(10**row_data['density_chi2']*u.cm**-3,
            #                                             width,
            #                                             lengthscale,
            #                                             width/lengthscale,
            #                                             columns['higaldusttem'][-1]*u.K,
            #                                             crir=0./u.s))

        if len(set(len(c) for k,c in columns.items())) != 1:
            print("Columns are different lengths.  This is not allowed.")
            import ipdb; ipdb.set_trace()

        for c in columns:
            assert len(columns[c]) == ii+1

        if plot_some and not is_bad and (ii-nbad % 100 == 0 or ii-nbad < 50):
            try:
                log.info("T: [{tmin1sig_chi2:7.2f},{expected_temperature:7.2f},{tmax1sig_chi2:7.2f}]"
                         "  R={ratio321303:8.4f}+/-{eratio321303:8.4f}"
                         "  Smean303={Smean303:8.4f} +/- {e303:8.4f}"
                         "  Stot303={Stot303:8.2e}  npix={npix:6d}"
                         .format(Smean303=Smean303, Stot303=Stot303,
                                 npix=npix, e303=error, **row_data))

                pl.figure(1)
                pl.clf()
                mf.denstemplot()
                pl.savefig(fpath("dendrotem/diagnostics/{0}_{1}.png".format(suffix,ii)))
                pl.figure(2).clf()
                mf.parplot1d_all(levels=[0.68268949213708585])
                pl.savefig(fpath("dendrotem/diagnostics/1dplot{0}_{1}.png".format(suffix,ii)))
                pl.draw()
                pl.show()
            except Exception as ex:
                print ex
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

    for mid,lo,hi,letter in (('expected_temperature','tmin1sig_chi2','tmax1sig_chi2','t'),
                             ('expected_density','dmin1sig_chi2','dmax1sig_chi2','d'),
                             ('expected_column','cmin1sig_chi2','cmax1sig_chi2','c')):
        catalog.add_column(table.Column(name='elo_'+letter,
                                        data=catalog[mid]-catalog[lo]))
        catalog.add_column(table.Column(name='ehi_'+letter,
                                        data=catalog[hi]-catalog[mid]))

    if write:
        catalog.write(tpath('PPV_H2CO_Temperature{0}.ipac'.format(suffix)), format='ascii.ipac')

    # Note that there are overlaps in the catalog, which means that ORDER MATTERS
    # in the above loop.  I haven't yet checked whether large scale overwrites
    # small or vice-versa; it may be that both views of the data are interesting.
    tcube = SpectralCube(data=tcubedata, wcs=cube303.wcs,
                         mask=cube303.mask, meta={'unit':'K'},
                         header=cube303.header,
                        )
    tcubeleaf = SpectralCube(data=tcubeleafdata, wcs=cube303.wcs,
                         mask=cube303.mask, meta={'unit':'K'},
                         header=cube303.header,
                        )

    if write:
        log.info("Writing TemperatureCube")
        outpath = 'TemperatureCube_DendrogramObjects{0}.fits'
        tcube.write(hpath(outpath.format(suffix)),
                    overwrite=True)

        outpath_leaf = 'TemperatureCube_DendrogramObjects{0}_leaves.fits'
        tcubeleaf.write(hpath(outpath_leaf.format(suffix)),
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

def do_321_dendro_temperatures_sharp():
    from dendrograms import dend321
    catalog, tcube = measure_dendrogram_properties(dend=dend321,
                                                   cube303=cube303,
                                                   cube321=cube321,
                                                   sncube=sncube, suffix="",
                                                   line='321',
                                                   write=False)
    catalog.write(tpath('PPV_H2CO_Temperature_321selected.ipac'), format='ascii.ipac')
    tcube.write(hpath('TemperatureCube_Dendrogram321Objects.fits'), overwrite=True)

    return catalog,tcube

def do_321_dendro_temperatures_smooth():
    from dendrograms import dend321sm
    assert sncubesm._wcs is cube303sm._wcs
    catalog, tcube = measure_dendrogram_properties(dend=dend321sm,
                                                   cube303=cube303sm,
                                                   cube321=cube321sm,
                                                   cube13co=cube13cosm,
                                                   cube18co=cube18cosm,
                                                   noise_cube=sm_noise_cube,
                                                   sncube=sncubesm,
                                                   suffix="_smooth",
                                                   line='321',
                                                   write=False)

    catalog.write(tpath('PPV_H2CO_Temperature_321selected_smooth.ipac'), format='ascii.ipac')
    tcube.write(hpath('TemperatureCube_Dendrogram321Objects_smooth.fits'), overwrite=True)

    return catalog, tcube

def do_321_dendro_temperatures_both():
    do_321_dendro_temperatures_sharp()
    do_321_dendro_temperatures_smooth()
