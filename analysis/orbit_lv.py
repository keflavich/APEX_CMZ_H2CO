import numpy as np
from spectral_cube import SpectralCube, BooleanArrayMask
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

pixscale = wcs.utils.proj_plane_pixel_scales(cube.wcs.celestial).mean()*u.deg
twoarcmin = int((2*u.arcmin / pixscale).decompose())
threearcmin = int((3*u.arcmin / pixscale).decompose())
fourarcmin = int((4*u.arcmin / pixscale).decompose())

def make_lv(minusheight=threearcmin, plusheight=threearcmin, ):

    print("minus = {0}, plus = {1}".format(minusheight, plusheight))

    mask = np.zeros(cube.shape[1:], dtype='bool')

    for lon, lat in zip(longitudes,latitudes):
        if np.isfinite(lat) and np.isfinite(lon):
            xpix, ypix = map(int, map(np.round, cube.wcs.celestial.wcs_world2pix(lon, lat, 0)))
            mask[ypix-minusheight:ypix+plusheight, xpix] = True

    print("mask sum = {0}".format(mask.sum()))

    mcube = cube.with_mask(BooleanArrayMask(mask,
                                            wcs=cube.wcs,
                                            shape=cube.shape)).minimal_subcube(spatial_only=True)
    lv = mcube.sum(axis=1)

    fig = pl.figure()
    w = lv.wcs
    w.wcs.ctype[0] = '_GLONCAR'
    ax = pl.subplot(projection=w)
    ax.imshow(lv, interpolation='none', origin='lower', cmap='viridis')
    ax.set_xlabel('Galactic Longitude [deg]')
    ax.set_ylabel('$V_{LSR}$ [km s$^{-1}$]')
    fig.savefig(fpath("lvdiagram_13CO_{0}_to_{1}_pixels.pdf".format(minusheight, plusheight)))


    mx = cube.max(axis=0)
    mx.quicklook()
    mx.FITSFigure.show_contour(mask.astype('int'), levels=[0.5], colors=['r'])
    mx.FITSFigure.recenter(0,0,0.8)
    mx.FITSFigure.savefig(fpath("lvdiagram_cutout_region_13CO_{0}_to_{1}_pixels.pdf".format(minusheight, plusheight)))

    return lv

if __name__=='__main__':

    pl.close('all') # close open figures
    print("two, three, four arcmin = {0}, {1}, {2} pixels".format(twoarcmin, threearcmin, fourarcmin))
    lv33 = make_lv(minusheight=threearcmin, plusheight=threearcmin)
    lv24 = make_lv(minusheight=twoarcmin, plusheight=fourarcmin)
    lv22 = make_lv(minusheight=twoarcmin, plusheight=twoarcmin)
    lv22.write(paths.adpath('lv/lvdiagram_13CO_{0}_to_{1}_arcmin.fits'.format(2,2)))
    lv33.write(paths.adpath('lv/lvdiagram_13CO_{0}_to_{1}_arcmin.fits'.format(3,3)))
    lv24.write(paths.adpath('lv/lvdiagram_13CO_{0}_to_{1}_arcmin.fits'.format(2,4)))
