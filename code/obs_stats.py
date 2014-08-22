from paths import rpath,mpath
from astropy.io import fits
import pyregion

noise_img = fits.open(mpath('APEX_H2CO_merge_high_sub_noise.fits'))[0]
nhits_img = fits.open(mpath('APEX_H2CO_merge_high_nhits.fits'))[0]

map_regions = pyregion.open(rpath('targets_8x8.reg'))
