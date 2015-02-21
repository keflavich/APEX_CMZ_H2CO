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
from astropy.wcs import WCS
from astropy import table
import paths
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

import shapely.geometry as geom
import dendrograms
reload(dendrograms)
from dendrograms import dend, dendsm, catalog, catalogsm

# obsolete x,y = np.loadtxt(apath('orbit_K14.dat')).T
table = ascii.read(apath('orbit_K14_2.dat'), format='basic', comment="#", guess=False) 
coords = coordinates.SkyCoord(table['l']*u.deg, table['b']*u.deg, frame='galactic')
P = pvextractor.Path(coords, width=300*u.arcsec)

dl = (table['l'][1:]-table['l'][:-1])
db = (table['b'][1:]-table['b'][:-1])
dist = (dl**2+db**2)**0.5
cdist = np.zeros(dist.size+1)
cdist[1:] = dist.cumsum()

vmin,vmax=10,150
npoints = 10000

offset = np.linspace(0, cdist.max(), npoints)
timeline = np.interp(offset, cdist, table['t'])
veloline = np.interp(offset, cdist, table["v'los"])

line = geom.LineString(table['l','b'])
full_line_points = [line.interpolate(x) for x in np.linspace(0, line.length, npoints)]
full_line_xy = np.array([[p.x, p.y] for p in full_line_points])

def nearest_point(l, b, v=None, points=full_line_xy, velos=veloline,
                  return_distance=False):

    x,y = points.T
    z = velos

    distance = ((l-x)**2 + (b-y)**2)**0.5
    if return_distance:
        return distance.min()
    
    return np.argmin(distance)



def offset_to_point(glon, glat):
    """
    Determine the offset along the orbit to the nearest point on an orbit to
    the specified point
    """

    line = geom.LineString(table['l','b'])
    if hasattr(glon, '__len__'):
        distances = [line.project(geom.Point(l,b))
                  for l,b in zip(glon,glat)]
        return distances
    else:
        point = geom.Point(glon, glat)
        return line.project(point)

#for cat,smooth in zip((catalog, catalogsm),("","_smooth")):
for cat,smooth in ((catalog,"",),):
    glon = cat['x_cen'] - (cat['x_cen'] > 180)*360
    glat = cat['y_cen']
    vcen = cat['v_cen']/1e3
    nearest_points = [nearest_point(l,b,v) for l,b,v in zip(glon,glat,vcen)]
    model_xy = np.array([full_line_xy[iinp,:] for iinp in nearest_points])
    time = np.array([timeline[iinp] for iinp in nearest_points])
    modelvelo = np.array([veloline[iinp] for iinp in nearest_points])
    temperature = cat['temperature_chi2']
    distance = np.array([nearest_point(l,b,return_distance=True) for l,b in zip(glon,glat)])

    reftime = -2
    bricktime = 0.3
    cat.add_column(table.Column(time-reftime, name='OrbitTime'))
    cat.add_column(table.Column(modelvelo, name='ModelVelo'))
    cat.add_column(table.Column(distance, name='DistanceFromOrbit'))

    is_leaf = (cat['is_leaf'] == 'True')
    selection = ((np.abs((modelvelo-vcen))<35) &
                 (distance < 0.15) &
                 (vcen > -80))

    # make sure the Brick is included
    assert dend.structure_at((195,142,730)).idx in catalog['_idx'][selection]

    fig3 = pl.figure(3, figsize=(14,8))
    fig3.clf()
    ax3 = fig3.gca()

    segmentdata = {'alpha': [(0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
                   'blue': [(0.0, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
                   'green': [(0.0, 0.0, 0.0), (0.5, 0.75, 0.75), (1.0, 0.0, 0.0)],
                   'red': [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 1.0, 1.0)],
                   'alpha':[(0.0,0.5,0.5), (0.5, 0.5, 0.5), (1.0, 0.5, 0.5)],
                  }
    cmap = matplotlib.colors.LinearSegmentedColormap(name='rgb',
                                                   segmentdata=segmentdata)

    cbvmin,cbvmax = -50, 120
    color = cmap((vcen-cbvmin)/(cbvmax-cbvmin))
    cblabel = r"$v_{LSR}$ (km s$^{-1}$)"

    sc = ax3.scatter(time[selection&is_leaf]-reftime,
                     temperature[selection&is_leaf], marker='.',
                     c=color[selection&is_leaf], #alpha=0.5,
                     s=5000*cat['Smean303'][selection&is_leaf], edgecolor='none')
    ax3.errorbar(time[selection&is_leaf]-reftime,
                 temperature[selection&is_leaf],
                 yerr=[cat['temperature_chi2'][selection&is_leaf]-cat['tmin1sig_chi2'][selection&is_leaf],
                       -(cat['temperature_chi2'][selection&is_leaf]-cat['tmax1sig_chi2'][selection&is_leaf]),],
                 marker=None, linestyle='none', color='k', alpha=0.1,
                 capsize=0)

    sc3 = ax3.scatter(timeline-reftime, 0*timeline+5,
                      c=cmap((veloline-cbvmin)/(cbvmax-cbvmin)),
                      marker='|', edgecolor=cmap((veloline-cbvmin)/(cbvmax-cbvmin)),
                      s=100, alpha=0.5,
                     )

    # plot the non-leaf objects too
    sc2 = ax3.scatter(time[selection&~is_leaf]-reftime,
                     temperature[selection&~is_leaf], marker='.',
                     c=color[selection&~is_leaf], alpha=0.2,
                     s=5000*cat['Smean303'][selection&~is_leaf],
                     edgecolor='none')
    ax3.set_xlabel("Time since 1$^\\mathrm{st}$ pericenter passage [Myr]", size=24, labelpad=10)

    ytext=180
    ax3.text(bricktime, ytext, "Brick", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='k', weight='bold')
    ax3.text(bricktime+0.43, ytext, "Sgr B2", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='k', weight='bold')
    pl.setp(ax3.get_xticklabels(), fontsize=20)
    pl.setp(ax3.get_yticklabels(), fontsize=20)
    ax3.set_ylim(0, 200)
    ax3.set_ylabel("Temperature (K)")

    sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=cbvmin, vmax=cbvmax),
                                      cmap=cmap)
    sm._A = []
    cb = fig3.colorbar(sm)
    cb.set_label(cblabel)

    if smooth:
        ax3.axis([-0.1,1,0,200])
    else:
        ax3.axis([-0.0,1,0,200])
    pl.savefig(fpath("orbits/dendro{0}_temperature_vs_time_firstMyr.pdf".format(smooth)),
               bbox_inches='tight')

    L, = ax3.plot([0,0.8], [30,100], 'k--', linewidth=3, alpha=0.3, zorder=-10)
    pl.savefig(fpath("orbits/dendro{0}_temperature_vs_time_firstMyr_trendline.pdf".format(smooth)),
               bbox_inches='tight')
    L.set_visible(False)

    if smooth:
        ax3.set_xlim(-0.1,4.6)
    else:
        ax3.set_xlim(-0.0,4.6)
    ax3.text(bricktime+3.58, ytext*135./140., "20 km s$^{-1}$", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='k', weight='bold')
    ax3.text(bricktime+3.66, ytext*135./140., "50 km s$^{-1}$", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='k', weight='bold')
    ax3.text(bricktime+3.28, ytext, "Sgr C", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='k', weight='bold')

    pl.savefig(fpath("orbits/dendro{0}_temperature_vs_time.pdf".format(smooth)),
               bbox_inches='tight')

    L.set_visible(True)
    pl.savefig(fpath("orbits/dendro{0}_temperature_vs_time_trendline.pdf".format(smooth)),
               bbox_inches='tight')

    fig4 = pl.figure(4, figsize=(14,8))
    fig4.clf()
    ax4 = fig4.gca()
    ax4.hist(cat['temperature_chi2'][selection&is_leaf], histtype='stepfilled', edgecolor='none', alpha=0.5)
    ax4.hist(cat['temperature_chi2'][(~selection)&is_leaf], histtype='stepfilled', edgecolor='none', alpha=0.5)
