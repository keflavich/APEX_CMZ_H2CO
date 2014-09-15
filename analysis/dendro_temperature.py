from spectral_cube import SpectralCube, BooleanArrayMask
import numpy as np
import FITS_tools
from numpy.lib.stride_tricks import as_strided
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy import log
from astropy.utils.console import ProgressBar
from astrodendro import Dendrogram

from paths import hpath,mpath
from constrain_parameters import paraH2COmodel
from masked_cubes import cube303m,cube321m
from noise import noise, noise_cube
from higal_gridded import column_regridded
from common_constants import logabundance,elogabundance

mf = paraH2COmodel()
dend = Dendrogram.load_from(hpath("DendroMask_H2CO303202_signal_to_noise.hdf5"))

biggest_tree = dend[89]

metadata = {}
metadata['data_unit'] = u.K
metadata['spatial_scale'] =  7.2 * u.arcsec
metadata['beam_major'] =  30 * u.arcsec
metadata['beam_minor'] =  30 * u.arcsec
metadata['wavelength'] =  218.22*u.GHz
catalog = ppv_catalog(biggest_tree.descendants, metadata)


noise_flat = as_strided(noise, shape=mask.shape, strides=(0,)+noise.shape)[mask]

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
        's303',
        's321',
        'e303',
        'e321',
        'r303321',
        'er303321',
]
columns = {k:[] for k in keys}

pb = ProgressBar(len(biggest_tree.descendants))
for ii,structure in enumerate(biggest_tree.descendents):
    mask = BooleanArrayMask(structure.get_mask(), wcs=cube303m.wcs)
    view = cube303m.subcube_slices_from_mask(mask)
    submask = mask[view]

    c303 = cube303m[view].with_mask(submask)
    c321 = cube303m[view].with_mask(submask)

    s303 = c303.sum()
    s321 = c321.sum()
    r321303 = s321/s303

    error = (noise_cube[view][submask]).sum() / mask.sum()**0.5
    var = error**2
    er303321 = (r321303**2 * (var/s303**2 + var/s321**2))**0.5

    columns['s303'].append(s303)
    columns['s321'].append(s321)
    columns['e303'].append(error)
    columns['e321'].append(error)
    columns['r303321'].append(r321303)
    columns['er303321'].append(er321303)

    mask2d = mask.max(axis=0)[view[1:]]
    h2column = np.log10(column_regridded[view[1:]][mask2d].mean() * 1e22)
    elogh2column = elogabundance

    mf.set_constraints(ratio303321=r321303, eratio303321=er321303,
                       #ratio321322=ratio2, eratio321322=eratio2,
                       logh2column=logh2column, elogh2column=elogh2column,
                       logabundance=logabundance, elogabundance=elogabundance,
                       taline303=s303, etaline303=error,
                       taline321=s321, etaline321=error,
                       linewidth=5)
    row_data = mf.get_parconstraints()

    for k in row_data:
        columns[k].append(row_data[k])

    if ii % 100 == 0 or ii < 50:
        log.info("T: [{tmin1sig_chi2:7.2f},{temperature_chi2:7.2f},{tmax1sig_chi2:7.2f}]  R={ratio303321:6.2f}+/-{eratio303321:6.2f}".format(**row_data))
    else:
        pb.update(ii+1)
