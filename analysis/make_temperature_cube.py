from spectral_cube import SpectralCube, BooleanArrayMask
import numpy as np
import FITS_tools
from numpy.lib.stride_tricks import as_strided
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy import log
from astropy.utils.console import ProgressBar

from paths import hpath
from constrain_parameters import paraH2COmodel

nsigma = 5 # start big to minimize # of failures

mf = paraH2COmodel()

cube303 = SpectralCube.read(hpath('APEX_H2CO_303_202.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
cube321 = SpectralCube.read(hpath('APEX_H2CO_321_220.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
mask = (fits.getdata(hpath('APEX_H2CO_303_202_mask.fits')).astype('bool') &
        cube303.mask.include(cube303._data, cube303.wcs) &
        cube321.mask.include(cube321._data, cube321.wcs))
bmask = BooleanArrayMask(mask, cube303.wcs)
cube303m = cube303.with_mask(bmask)
cube321m = cube321.with_mask(bmask)

column_image = fits.open('/Users/adam/work/gc/gcmosaic_column_conv36.fits')[0]
dusttem_image = fits.open('/Users/adam/work/gc/gcmosaic_temp_conv36.fits')[0]

apex_header = cube303[0,:,:].hdu.header
column_regridded = FITS_tools.hcongrid.hcongrid_hdu(column_image, apex_header)

noise = fits.getdata(hpath('APEX_H2CO_303_202_noise.fits'))
noise = cube303[:50].std(axis=0).value
noise_flat = as_strided(noise, shape=mask.shape, strides=(0,)+noise.shape)[mask]
var_flat = noise_flat**2

ratio303321 = (cube303m.flattened() / cube321m.flattened()).value
eratio303321 = (ratio303321**2 * (var_flat/cube303m.flattened().value**2 + var_flat/cube321m.flattened().value**2))**0.5

indices = np.where(mask)
usable = (eratio303321*nsigma < ratio303321) & (eratio303321 > 0) & (ratio303321 > 0) & (noise_flat > 1e-10) & (noise_flat < 10) & (ratio < 100)
ngood = np.count_nonzero(usable)
usable_indices = [ind[usable] for ind in indices]
uz,uy,ux = usable_indices

column_flat = column_regridded.data[uy,ux]
uratio303321 = ratio303321[usable]
ueratio303321 = eratio303321[usable]
utline303 = cube303.flattened()[usable]
utline321 = cube321.flattened()[usable]
unoise = noise_flat[usable]

log.info("Out of {size} values, {ngood} are usable "
         "for fitting ({pct:0.1f}%).".format(size=ratio303321.size,
                                             ngood=ngood,
                                             pct=100*ngood/float(ratio303321.size)))

tcube = np.empty_like(mask, dtype='float')
tcube[:] = np.nan

logabundance = np.log10(1.2e-9)
elogabundance = 1.0
elogh2column = elogabundance
linewidth = 5 # this is ugly...

pb = ProgressBar(ngood)
for ii,((z,y,x),rat,erat,col,ta303,ta321,err) in enumerate(zip(zip(*usable_indices),
                                                           uratio303321, ueratio303321,
                                                           column_flat, utline303,
                                                           utline321, unoise)):
    logh2column = np.log10(col)+22

    mf.set_constraints(ratio303321=rat, eratio303321=erat,
                       #ratio321322=ratio2, eratio321322=eratio2,
                       logh2column=logh2column, elogh2column=elogh2column,
                       logabundance=logabundance, elogabundance=elogabundance,
                       taline303=ta303.value, etaline303=err,
                       taline321=ta321.value, etaline321=err,
                       linewidth=linewidth)
    row_data = mf.get_parconstraints()
    tcube[z,y,x] = row_data['temperature_chi2']
    row_data['ratio303321'] = rat
    row_data['eratio303321'] = erat

    if ii % 100 == 0 or ii < 50:
        log.info("T: [{tmin1sig_chi2:7.2f},{temperature_chi2:7.2f},{tmax1sig_chi2:7.2f}]  R={ratio303321:6.2f}+/-{eratio303321:6.2f}".format(**row_data))
    else:
        pb.update(ii)

print()
