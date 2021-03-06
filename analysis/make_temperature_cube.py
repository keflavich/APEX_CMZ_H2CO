import os
from spectral_cube import SpectralCube, BooleanArrayMask
import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy import log
from astropy.utils.console import ProgressBar

from paths import hpath
from constrain_parameters import paraH2COmodel
from masked_cubes import cube303m,cube321m,cube303,cube321,mask
from noise import noise, noise_cube
from common_constants import logabundance,elogabundance
from higal_gridded import column_regridded
from ratio_cubes import ratio303321, eratio303321, noise_flat

log.warn("Extremely slow for a very small fraction of overall data."
         "Try a different method: make_piecewise_temcube.py")

nsigma = 5 # start big to minimize # of failures

mf = paraH2COmodel()

indices = np.where(mask)
usable = ((eratio303321*nsigma < ratio303321) & (eratio303321 > 0) &
          (ratio303321 > 0) & (noise_flat > 1e-10) & (noise_flat < 10) &
          (ratio303321 < 100))
assert len(usable) == len(indices)
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

log.info("Using temperature memmap.")
if os.path.exists('temperature_memmap'):
    mode = 'r+'
else:
    mode = 'w+'
tcube = np.memmap('temperature_memmap', mode=mode, dtype='float32', shape=mask.shape)
nfilled, nbad = np.count_nonzero(tcube), np.count_nonzero(np.isnan(tcube))
log.info("There are {0} already-filled entries and {1} nan"
         " entries.".format(nfilled,nbad))
if nfilled + nbad == tcube.size:
    log.warn("The temperature memmap is already full.  To recreate it,"
             " delete the file then re-run this code.")

elogh2column = elogabundance
linewidth = 5 # this is ugly...

pb = ProgressBar(ngood)
for ii,((z,y,x),rat,erat,col,ta303,ta321,err) in enumerate(zip(zip(*usable_indices),
                                                               uratio303321, ueratio303321,
                                                               column_flat, utline303,
                                                               utline321, unoise)):
    if tcube[z,y,x] == 0:
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
        tcube.flush()
    else:
        pb.update(ii)

tcube[tcube==0] = np.nan
tCube = SpectralCube(tcube, cube303.wcs, mask=BooleanArrayMask(np.isfinite(tcube), wcs=cube303.wcs))
tCube.write(hpath('chi2_temperature_cube.fits'), overwrite=True)

print()
