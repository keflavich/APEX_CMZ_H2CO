import pyspeckit
import numpy as np
from astropy.io import fits
from pyspeckit.spectrum import models
from pyspeckit.spectrum.models.model import SpectralModel
import FITS_tools

# create the Formaldehyde Radex fitter
# This step cannot be easily generalized: the user needs to read in their own grids
texgrid1 = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321-220_5kms_temperature_para_tex1.fits')
taugrid1 = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321-220_5kms_temperature_para_tau1.fits')
texgrid2 = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321-220_5kms_temperature_para_tex2.fits')
taugrid2 = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321-220_5kms_temperature_para_tau2.fits')
hdr = fits.getheader('/Users/adam/work/h2co/radex/thermom/303-202_321-220_5kms_temperature_para_tau2.fits')

texgrid1b = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_5kms_temperature_para_tex1.fits')
taugrid1b = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_5kms_temperature_para_tau1.fits')
texgrid2b = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_5kms_temperature_para_tex2.fits')
taugrid2b = fits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_5kms_temperature_para_tau2.fits')
hdrb = fits.getheader('/Users/adam/work/h2co/radex/thermom/303-202_322-221_5kms_temperature_para_tau2.fits')

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
                                  parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)],
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
    cube.Registry.add_fitter('h2co_simple', simple_fitter, 4, multisingle='multi')

    import pyregion

    regs = pyregion.open(regpath+'spectral_apertures.reg')

    for reg in regs:
        name = reg.attr[1]['text']
        if name not in spectra:
            sp = cube.get_apspec(reg.coord_list,coordsys='galactic',wunit='degree')
            sp.specname = reg.attr[1]['text']
            sp.error[:] = sp.stats((218e9,218.1e9))['std']
        sp.plotter()
        #sp.specfit(fittype='h2co_mm_radex', multifit=True, guesses=[100,14.3,4.0,40,7.0],
        #           limits=[(20,200),(11,15),(3,5.5),(-105,105),(2,18)], limited=[(True,True)]*5,
        #           fixed=[False,True,True,False,False], quiet=False,)
        sp.specfit(fittype='h2co_simple', multifit=True,
                   guesses=[1,25,5,0.5,1])
        sp.specfit.plot_fit()
        spectra[sp.specname] = sp

    cube.fiteach(fittype='h2co_simple', guesses=[1,25,5,0.5,1], multicore=8, errmap=noise, sigmacut=5)
    amp,vel,wid,ratio,ch3oh = cube.parcube
    eamp,evel,ewid,eratio,ech3oh = cube.errcube
    ok = (amp > 0) & (amp > eamp*2) & (vel > -100) & (vel < 150) & (wid > ewid*2) & (ratio > 0) & (ratio < 1) & (eratio < 0.3)
    hdr = FITS_tools.strip_headers.flatten_header(cube.header)
    ratio[True-ok] = np.nan
    ratioimg = fits.PrimaryHDU(data=ratio, header=hdr)
    ratioimg.writeto('H2CO_fitted_ratios.fits', clobber=True)
    from shfi_otf_pipeline.make_apex_cubes import tm
    tmap = tm(ratio)
    ratioimg.data = tmap
    ratioimg.writeto('H2CO_fitted_tmap.fits', clobber=True)
    # use_nearest_as_guess=True, 


    individual_fits=False
    if individual_fits:
        flux = 3.6 # Jy
        col_per_jy = 2e22 # cm^-2
        dvdpc = 5.0 # km/s/pc
        logX = -8.3
        logcol = np.log10(flux*col_per_jy/dvdpc) + logX

        spectra['WarmSpot'].specfit(fittype='h2co_mm_radex', multifit=True,
                                    guesses=[100,logcol,4.5,35,3.0],
                                    limits=[(20,200),(11,15),(3,5.5),(-105,105),(1,5)],
                                    limited=[(True,True)]*5,
                                    fixed=[False,True,True,False,False],
                                    quiet=False,)

        spectra['WarmSpot'].specfit(fittype='h2co_mm_radex', multifit=True,
                                    guesses=[100,logcol+np.log10(2/3.),4.5,27,3.0]+[100,logcol+np.log10(1.0/3.),4.5,53,3.0],
                                    limits=[(20,200),(11,15),(3,5.5),(-105,105),(1,8)]+[(20,200),(11,15),(3,5.5),(-105,105),(1,6)],
                                    limited=[(True,True)]*10,
                                    fixed=[False,True,True,False,False]*2,
                                    quiet=False,)

        flux = 5.0
        logX = -8.5
        logcol = np.log10(flux*col_per_jy/dvdpc / 2.) + logX

        spectra['Brick SW'].specfit(fittype='h2co_mm_radex', multifit=True,
                                    guesses=[133,12.94,5.977,37.17,9.88],
                                    limits=[(20,200),(11,15),(3,6.5),(-105,105),(1,15)],
                                    limited=[(True,True)]*5,
                                    #fixed=[False,False,True,False,False],
                                    fixed=[False,True,True,False,False])

        spectra['50kmsColdExtension'].specfit(fittype='h2co_mm_radex', multifit=True,
                                    guesses=[33,12.94,4.0,24.4,7],
                                    limits=[(20,200),(11,15),(3,6.5),(-105,105),(1,15)],
                                    limited=[(True,True)]*5,
                                    #fixed=[False,False,True,False,False],
                                    fixed=[False,False,True,False,False])

        #spectra['Sgr B2 SW'].specfit(fittype='h2co_mm_radex', multifit=True,
        #                            guesses=[125,14.14,4.0,47.72,5.66]+[125,14.14,4.0,55.72,5.66],
        #                            limits=[(20,200),(11,15),(3,6.5),(-105,105),(1,15)]*2,
        #                            limited=[(True,True)]*5*2,
        #                            #fixed=[False,False,True,False,False],
        #                            fixed=[False,False,True,True,True]+[False,False,True,False,False])



    dopymc = False
    if dopymc:
        import agpy

        sp = spectra['20 kms']
        # SHOULD BE 
        spmc = sp.specfit.get_pymc(use_fitted_values=True, db='hdf5', dbname='h2co_mm_fit_20kmsCld.hdf5')
        #spmc = sp.specfit.fitter.get_pymc(sp.xarr, sp.data, sp.error,
        #                                  use_fitted_values=True, db='hdf5',
        #                                  dbname='h2co_mm_fit_20kmsCld.hdf5')
        spmc.sample(100000)
        agpy.pymc_plotting.hist2d(spmc, 'TEMPERATURE0', 'DENSITY0', doerrellipse=False, clear=True, bins=50, fignum=4)
        agpy.pymc_plotting.hist2d(spmc, 'TEMPERATURE0', 'COLUMN0', doerrellipse=False, clear=True, bins=50,fignum=5)
        agpy.pymc_plotting.hist2d(spmc, 'DENSITY0', 'COLUMN0', doerrellipse=False, clear=True, bins=50,fignum=6)


        pars = dict([(k,spmc.trace(k)[-50:]) for k in sp.specfit.parinfo.keys()])
        sp.plotter.autorefresh=False
        for ii in xrange(0,50):
            sp.specfit.plot_model([pars[k][ii] for k in sp.specfit.parinfo.keys()],
                                  clear=False,
                                  composite_fit_color='r',
                                  plotkwargs={'alpha':0.01})

        sp.plotter.refresh()
