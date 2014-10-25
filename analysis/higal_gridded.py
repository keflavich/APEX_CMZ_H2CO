import FITS_tools
from astropy.io import fits
from masked_cubes import cube303
from astropy import convolution

column_image = fits.open('/Users/adam/work/gc/gcmosaic_column_conv36.fits')[0]
dusttem_image = fits.open('/Users/adam/work/gc/gcmosaic_temp_conv36.fits')[0]

# fix NaNs by convolving
col_conv = convolution.convolve_fft(column_image.data,
                                    convolution.Gaussian2DKernel(2),
                                    interpolate_nan=True)
whnan = np.isnan(column_image.data)
column_image.data[whnan] = col_conv[whnan]

dusttem_conv = convolution.convolve_fft(dusttem_image.data,
                                    convolution.Gaussian2DKernel(2),
                                    interpolate_nan=True)
whnan = np.isnan(dusttem_image.data)
dusttem_image.data[whnan] = dusttem_conv[whnan]


apex_header = cube303[0,:,:].hdu.header
column_regridded = FITS_tools.hcongrid.hcongrid_hdu(column_image, apex_header)
