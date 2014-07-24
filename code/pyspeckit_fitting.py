import pyspeckit
import numpy as np
from astropy.io import fits
from pyspeckit.spectrum import models
from pyspeckit.spectrum.models.model import SpectralModel
import FITS_tools
from paths import h2copath, mergepath, figurepath
from shfi_otf_pipeline.make_apex_cubes import tm

# create the Formaldehyde Radex fitter
# This step cannot be easily generalized: the user needs to read in their own grids
texgrid1 = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321-220_temperature_para_300Kmax_5kms_tex1.fits')
taugrid1 = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321-220_temperature_para_300Kmax_5kms_tau1.fits')
texgrid2 = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321-220_temperature_para_300Kmax_5kms_tex2.fits')
taugrid2 = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321-220_temperature_para_300Kmax_5kms_tau2.fits')
hdr = fits.getheader('/Users/adam/work/h2co/radex/thermom/303-202_321-220_temperature_para_300Kmax_5kms_tau2.fits')

texgrid1b = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_300Kmax_5kms_tex1.fits')
taugrid1b = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_300Kmax_5kms_tau1.fits')
texgrid2b = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_300Kmax_5kms_tex2.fits')
taugrid2b = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_300Kmax_5kms_tau2.fits')
hdrb = fits.getheader('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_300Kmax_5kms_tau2.fits')

# # this deserves a lot of explanation:
# # models.formaldehyde.formaldehyde_radex is the MODEL that we are going to fit
# # models.model.SpectralModel is a wrapper to deal with parinfo, multiple peaks,
# # and annotations
# # all of the parameters after the first are passed to the model function 

h2co_radex_fitter = SpectralModel(models.formaldehyde_mm.formaldehyde_mm_radex,
                                  5,
                                  parnames=['temperature','column','density','center','width'],
                                  parvalues=[50,12,4.5,0,1],
                                  parlimited=[(True,True), (True,True),
                                              (True,True), (False,False),
                                              (True,False)],
                                  parlimits=[(5,300), (11,17), (3,7), (0,0), (0,0)],
                                  parsteps=[0.01,0.01,0.1,0,0], fitunits='Hz',
                                  texgrid=((218.1,218.3,texgrid1b),
                                           (218.35,218.55,texgrid2b),
                                           (218.6,218.85,texgrid2)),
                                  # specify the frequency range over which the grid is valid (in GHz)
                                  taugrid=((218.1,218.3,taugrid1b),
                                           (218.35,218.55,taugrid2b),
                                           (218.6,218.85,taugrid2)),
                                  hdr=hdrb,
                                  shortvarnames=("T","N","n","v","\\sigma"),
                                  # specify the parameter names (TeX is OK)
                                  grid_vwidth=5.0,
                                  )


def simplemodel(xarr, amplitude, velocity, width, ratio, amp2):

    x = xarr.as_unit('km/s')
    voff1 = x.x_to_coord(218.475632e9,'Hz')
    voff2 = x.x_to_coord(218.760066e9,'Hz')
    voff3 = x.x_to_coord(218.44005,'GHz')
    x = np.array(x) # make sure xarr is no longer a spectroscopic axis
    G = (amplitude*(np.exp(-(x-velocity)**2/(2.0*width**2)) +
                    ratio*np.exp(-(x-velocity-voff1)**2/(2.0*width**2)) +
                    ratio*np.exp(-(x-velocity-voff2)**2/(2.0*width**2))) +
         amp2 * np.exp(-(x-velocity-voff3)**2/(2.0*width**2)))

    return G

simple_fitter = SpectralModel(simplemodel, 5,
                              parnames=['Amplitude','Velocity','Width','Ratio','AmpCH3OH'],
                              parvalues=[1,0,1,0.5,0.5],
                              parlimited=[(True,False),
                                          (False,False),
                                          (True,False),
                                          (True,True),
                                          (True,False)],
                              parlimits=[(0,0),(0,0),(0,0),(0,1),(0,0)],
                              shortvarnames=('A','v',r'\sigma','R','A_2'),
                              )



if __name__ == "__main__":
    from paths import mergepath,regpath

    if 'cube' not in locals():
        cube = pyspeckit.Cube(mergepath+'APEX_H2CO_merge_high_vsmoothds.fits')
        etamb = 0.75 # http://www.apex-telescope.org/telescope/efficiency/
        cube.cube /= etamb
        noise = fits.getdata(mergepath+'APEX_H2CO_merge_high_vsmoothds_noise.fits') / etamb
        spectra = {}

    cube.Registry.add_fitter('h2co_mm_radex', h2co_radex_fitter, 5,
                             multisingle='multi')
    cube.Registry.add_fitter('h2co_simple', simple_fitter, 5, multisingle='multi')

    cube.fiteach(fittype='h2co_simple', guesses=[1,25,5,0.5,1], multicore=8, errmap=noise, sigmacut=4)

    hdr = FITS_tools.strip_headers.flatten_header(cube.header)

    amp,vel,wid,ratio,ch3oh = cube.parcube
    eamp,evel,ewid,eratio,ech3oh = cube.errcube
    ratioimg = fits.PrimaryHDU(data=ratio, header=hdr)
    ratioimg.writeto('H2CO_fitted_ratios_raw.fits', clobber=True)

    ok = (amp > 0) & (amp > eamp*5) & (vel > -100) & (vel < 150) & (wid > ewid*2) & (ratio > 0) & (ratio < 1) & (eratio < 0.3) & (eratio < ratio) & (wid < 20)
    ratio = ratio.copy()
    ratio[True-ok] = np.nan
    ratioimg = fits.PrimaryHDU(data=ratio, header=hdr)
    ratioimg.writeto('H2CO_fitted_ratios.fits', clobber=True)
    tmap = tm(ratio)
    ratioimg.data = tmap
    ratioimg.writeto('H2CO_fitted_tmap.fits', clobber=True)
