import pyspeckit
import numpy as np
from astropy import log
from astropy.io import fits
from pyspeckit.spectrum import models
from pyspeckit.spectrum.models.model import SpectralModel
import FITS_tools
from paths import h2copath, mergepath, figurepath, gpath
from make_ratiotem_cubesims import tm
import os

try:
    # create the Formaldehyde Radex fitter
    # This step cannot be easily generalized: the user needs to read in their own grids
    texgrid303 = fits.getdata(gpath('fjdu_pH2CO_303_tex_5kms.fits'))
    taugrid303 = fits.getdata(gpath('fjdu_pH2CO_303_tau_5kms.fits'))
    texgrid321 = fits.getdata(gpath('fjdu_pH2CO_321_tex_5kms.fits'))
    taugrid321 = fits.getdata(gpath('fjdu_pH2CO_321_tau_5kms.fits'))
    texgrid322 = fits.getdata(gpath('fjdu_pH2CO_322_tex_5kms.fits'))
    taugrid322 = fits.getdata(gpath('fjdu_pH2CO_322_tau_5kms.fits'))
    hdr = hdrb = fits.getheader(gpath('fjdu_pH2CO_303_tex_5kms.fits'))

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
                                      texgrid=((218.1,218.3,texgrid303),
                                               (218.35,218.55,texgrid322),
                                               (218.6,218.85,texgrid321)),
                                      # specify the frequency range over which the grid is valid (in GHz)
                                      taugrid=((218.1,218.3,taugrid303),
                                               (218.35,218.55,taugrid322),
                                               (218.6,218.85,taugrid321)),
                                      hdr=hdrb,
                                      shortvarnames=("T","N","n","v","\\sigma"),
                                      # specify the parameter names (TeX is OK)
                                      grid_vwidth=5.0,
                                      )
except IOError as ex:
    log.exception("Could not read files from disk: cannot load H2CO RADEX fitter")
    log.exception(ex)


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

def simplemodel2(xarr, amplitude, velocity, width, ratio321303, ratio322321, amp2):

    x = xarr.as_unit('km/s')
    voff322 = x.x_to_coord(218.475632e9,'Hz')
    voff321 = x.x_to_coord(218.760066e9,'Hz')
    voff3 = x.x_to_coord(218.44005,'GHz')
    x = np.array(x) # make sure xarr is no longer a spectroscopic axis
    G = (amplitude*(np.exp(-(x-velocity)**2/(2.0*width**2)) +
                    np.exp(-(x-velocity-voff322)**2/(2.0*width**2))*ratio321303*ratio322321 +
                    np.exp(-(x-velocity-voff321)**2/(2.0*width**2))*ratio321303) +
         amp2 * np.exp(-(x-velocity-voff3)**2/(2.0*width**2)))

    return G

simple_fitter2 = SpectralModel(simplemodel2, 6,
                              parnames=['Amplitude', 'Velocity', 'Width',
                                        'Ratio321303x', 'Ratio322321x',
                                        'AmpCH3OH'], 
                              parvalues=[1,0,1,0.5,1.0,0.5],
                              parlimited=[(True,False),
                                          (False,False),
                                          (True,False),
                                          (True,True),
                                          (True,True),
                                          (True,False)],
                              parlimits=[(0,0),(0,0),(0,0),(0,1),(0.3,1.1),(0,0)],
                              shortvarnames=('A','v',r'\sigma','R1','R2','A_2'),
                              )



if __name__ == "__main__":
    from paths import mergepath,regpath

    if 'cube' not in locals():
        cube = pyspeckit.Cube(mergepath+'APEX_H2CO_merge_high_vsmoothds.fits')
        etamb = 0.67 # http://www.apex-telescope.org/telescope/efficiency/ then from class's ruze from Katharina's #'s
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
