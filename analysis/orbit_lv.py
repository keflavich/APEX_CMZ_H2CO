import numpy as np
from spectral_cube import SpectralCube
import pylab as pl
from paths import mpath,apath,fpath,molpath,hpath
from astropy import units as u
from astropy import coordinates
from astropy.io import ascii
from astropy import wcs
import paths
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

# obsolete x,y = np.loadtxt(apath('orbit_K14.dat')).T
table = ascii.read(apath('orbit_K14_2.dat'), format='basic', comment="#", guess=False)
coords = coordinates.SkyCoord(table['l']*u.deg, table['b']*u.deg, frame='galactic')

molecule = '13CO_2014_merge'
cube = SpectralCube.read(molpath('APEX_{0}.fits'.format(molecule)))

longitudes = cube.world[0,0,:][2]
longitudes[longitudes > 180*u.deg] -= 360*u.deg

first_switch = np.where(np.diff(coords.l.wrap_at(180*u.deg)) < 0)[0][0]

latitudes = np.interp(longitudes, coords.l.wrap_at(180*u.deg)[:first_switch+1],
                      coords.b[:first_switch+1], left=np.nan, right=np.nan)

mask = np.zeros(cube.shape[1:], dtype='bool')

pixscale = wcs.utils.proj_plane_pixel_scales(cube.wcs.celestial).mean()*u.deg
height = int((3*u.arcmin / pixscale).decompose())

for lon, lat in zip(longitudes,latitudes):
    if np.isfinite(lat) and np.isfinite(lon):
        xpix, ypix = map(int, map(np.round, cube.wcs.celestial.wcs_world2pix(lon, lat, 0)))
        mask[ypix-height:ypix+height, xpix] = True

mcube = cube.with_mask(mask).minimal_subcube(spatial_only=True)
lv = mcube.sum(axis=1)

lv.quicklook()
