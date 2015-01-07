import numpy as np
import pvextractor
import spectral_cube
import aplpy
import pylab as pl
import matplotlib
import copy
from paths import mpath,apath,fpath,molpath,hpath
from astropy import units as u
from astropy import coordinates
from astropy.io import ascii
from astropy import log

# obsolete x,y = np.loadtxt(apath('orbit_K14.dat')).T
table = ascii.read(apath('orbit_K14_2.dat'), format='basic', comment="#", guess=False) 
coords = coordinates.SkyCoord(table['l']*u.deg, table['b']*u.deg, frame='galactic')
P = pvextractor.Path(coords, width=120*u.arcsec)

dl = (table['l'][1:]-table['l'][:-1])
db = (table['b'][1:]-table['b'][:-1])
dist = (dl**2+db**2)**0.5
cdist = np.zeros(dist.size+1)
cdist[1:] = dist.cumsum()


molecules = ('13CO_2014_merge', 'C18O_2014_merge', 'H2CO_303_202_bl', 'SiO_54_bl')

filenames = [molpath('APEX_{0}.fits'.format(molecule))
            for molecule in molecules]

molecules = molecules + ('H2CO_DendrogramTemperature',
                         'H2CO_DendrogramTemperature_smooth')
filenames.append(hpath('TemperatureCube_DendrogramObjects_Piecewise.fits'))
filenames.append(hpath('TemperatureCube_DendrogramObjects_smooth_Piecewise.fits'))

cmap = copy.copy(pl.cm.RdYlBu_r)
cmap.set_bad((0.5,)*3)
cmap.set_under((0.9,0.9,0.9,0.5))

for molecule,fn in zip(molecules[-2:],filenames[-2:]):
    log.info(molecule)
    cube = spectral_cube.SpectralCube.read(fn)

    pv = pvextractor.extract_pv_slice(cube, P, respect_nan=False)

    fig1 = pl.figure(1, figsize=(14,8))
    fig1.clf()
    F = aplpy.FITSFigure(pv, figure=fig1)
    if 'Temperature' in fn:
        F.show_colorscale(cmap=cmap, aspect=5, vmin=20, vmax=200)
        # This is where it fails...
        #F.add_colorbar()
        #divider = make_axes_locatable(F._ax1)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #pl.colorbar(F._ax1.images[0], cax=cax)
    else:
        F.show_grayscale()
    F.show_lines(np.array([[cdist, table["v'los"]*1e3]]), zorder=1000, color='r',
                 linewidth=3, alpha=0.25)
    F.save(fpath('orbits/KDL2014_orbit_on_{0}.pdf'.format(molecule)))

    fig2 = pl.figure(2)
    pl.clf()
    F2 = aplpy.FITSFigure(cube.moment0().hdu, convention='calabretta', figure=fig2)
    if 'Temperature' in fn:
        F2.show_colorscale(cmap=cmap)
        F2.add_colorbar()
    else:
        F2.show_grayscale()

    patches = P.to_patches(1, ec='red', fc='none',
                           #transform=ax.transData,
                           clip_on=True, #clip_box=ax.bbox,
                           wcs=cube.wcs)
    for patch in patches:
        patch.set_linewidth(0.5)
        patch.set_alpha(0.5)
        patch.zorder = 50

    patchcoll = matplotlib.collections.PatchCollection(patches, match_original=True)
    patchcoll.zorder=10

    c = F2._ax1.add_collection(patchcoll)

    F2._rectangle_counter += 1
    rectangle_set_name = 'rectangle_set_' + str(F2._rectangle_counter)

    F2._layers[rectangle_set_name] = c
    F2.recenter(0, -0.03, width=1.5, height=0.3)

    F2.save(fpath('orbits/KDL2014_orbitpath_on_{0}.pdf'.format(molecule)))
