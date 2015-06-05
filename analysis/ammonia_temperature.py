from __future__ import print_function
from astropy import units as u
import numpy as np
import pyspeckit
from spectral_cube import SpectralCube
import paths
import pyregion

cube11 = SpectralCube.read(paths.adpath('G357.3-003.9-NH3-11-cube.fits'))
cube22 = SpectralCube.read(paths.adpath('G357.3-003.9-NH3-22-cube.fits'))
regions = pyregion.open(paths.rpath('target_fields_8x8.reg'))

sc11 = cube11.subcube_from_ds9region(regions)
sc22 = cube22.subcube_from_ds9region(regions)
sp_s11 = sc11.mean(axis=(1,2))
sp_s22 = sc22.mean(axis=(1,2))
print("Integrated line ratio 1-1/2-2: {0}".format(sp_s11.sum()/sp_s22.sum()))

filling_factor = 0.1

sp11 = pyspeckit.Spectrum(data=sp_s11.value/filling_factor, xarr=cube11.spectral_axis,
                          header=cube11.header,
                          xarrkwargs={'refX': cube11.wcs.wcs.restfrq*u.Hz,
                                      'velocity_convention': 'radio'})
sp22 = pyspeckit.Spectrum(data=sp_s22.value/filling_factor, xarr=cube22.spectral_axis,
                          header=cube22.header,
                          xarrkwargs={'refX': cube22.wcs.wcs.restfrq*u.Hz,
                                      'velocity_convention': 'radio'})

input_dict = {'oneone': sp11, 'twotwo': sp22}

T,F = True,False
# WHY are these out oflimits??!?!?!
raise Exception("There seems to be a bug in this next step; in any case it doesn't work.")
spf,spectra = pyspeckit.wrappers.fitnh3.fitnh3tkin(input_dict, dobaseline=False,
                                                   guesses=[55, 30, 12.5, 10,
                                                            50., 0.5],
                                                   fixed=[F,F,F,F,F,T],
                                                   limitedmin=[T]*6,
                                                   limitedmax=[T]*6,
                                                   limits=[(10,100),
                                                           (0,100),
                                                           (12, 16),
                                                           (3,30),
                                                           (45,55),
                                                           (0,1)],
                                                  )
pyspeckit.wrappers.fitnh3.plot_nh3(spf, spectra, show_hyperfine_components=True)
