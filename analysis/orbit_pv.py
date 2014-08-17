import numpy as np
import pvextractor
import spectral_cube
from astropy import units as u
from astropy import coordinates

x,y = np.loadtxt('orbit_K14.dat').T
coords = coordinates.SkyCoord(x*u.deg, y*u.deg, frame='galactic')
P = pvextractor.Path(coords)
cube = spectral_cube.SpectralCube.read('APEX_13CO_2014_merge.fits')

pv = pvextractor.extract_pv_slice(cube, P)

