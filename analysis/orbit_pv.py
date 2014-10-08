import numpy as np
import pvextractor
import spectral_cube
import aplpy
import pylab as pl
from paths import mpath,apath,fpath
from astropy import units as u
from astropy import coordinates
from astropy.io import ascii

x,y = np.loadtxt(apath('orbit_K14.dat')).T
table = ascii.read(apath('orbit_K14_2.dat'), format='basic', comment="#", guess=False) 
coords = coordinates.SkyCoord(table['l']*u.deg, table['b']*u.deg, frame='galactic')
P = pvextractor.Path(coords)
cube = spectral_cube.SpectralCube.read(mpath('APEX_13CO_2014_merge.fits'))

pv = pvextractor.extract_pv_slice(cube, P)

dl = (table['l'][1:]-table['l'][:-1])
db = (table['b'][1:]-table['b'][:-1])
dist = (dl**2+db**2)**0.5
cdist = np.zeros(dist.size+1)
cdist[1:] = dist.cumsum()

fig1 = pl.figure(1, figsize=(14,8))
fig1.clf()
F = aplpy.FITSFigure(pv, figure=fig1)
F.show_grayscale()
F.show_lines(np.array([[cdist, table["v'los"]*1e3]]), zorder=1000, color='r',
             linewidth=3, alpha=0.5)
F.save(fpath('KDL2014_orbit_on_APEX13CO2-1.pdf'))
