import FITS_tools
from astropy.io import fits
from masked_cubes import cube303

column_image = fits.open('/Users/adam/work/gc/gcmosaic_column_conv36.fits')[0]
dusttem_image = fits.open('/Users/adam/work/gc/gcmosaic_temp_conv36.fits')[0]

apex_header = cube303[0,:,:].hdu.header
column_regridded = FITS_tools.hcongrid.hcongrid_hdu(column_image, apex_header)
