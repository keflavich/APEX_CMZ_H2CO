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
import paths
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

# obsolete x,y = np.loadtxt(apath('orbit_K14.dat')).T
table = ascii.read(apath('orbit_K14_2.dat'), format='basic', comment="#", guess=False) 
coords = coordinates.SkyCoord(table['l']*u.deg, table['b']*u.deg, frame='galactic')
P = pvextractor.Path(coords, width=300*u.arcsec)

dl = (table['l'][1:]-table['l'][:-1])
db = (table['b'][1:]-table['b'][:-1])
dist = (dl**2+db**2)**0.5
cdist = np.zeros(dist.size+1)
cdist[1:] = dist.cumsum()


def offset_to_point(glon, glat):
    """
    Determine the offset along the orbit to the nearest point on an orbit to
    the specified point
    """
    import shapely.geometry as geom

    line = geom.LineString(table['l','b'])
    point = geom.Point(glon, glat)
    return line.project(point)
    point_on_line = line.interpolate(line.project(point))
    cl, cb = point_on_line.x, point_on_line.y

    # distance from point to each vertex
    vertex_distances = ((cl-table['l'])**2+(cb-table['b'])**2)
    closest_vertex = vertex_distances.argmin()
    after = vertex_distances[closest_vertex-1] < vertex_distances[closest_vertex+1]
    if after:
        i1 = closest_vertex
    else:
        i1 = closest_vertex-1
    i2 = i1+1

    line = geom.LineString(table['l','b'][i1:i2+1])
    point = geom.Point(glon, glat)
    point_on_line = line.interpolate(line.project(point))

    total_offset = ((((table['l'][1:i2] - table['l'][0:i1])**2 + 
                      (table['b'][1:i2] - table['b'][0:i1])**2)**0.5).sum() + 
                    ((point_on_line.x-table['l'][i2])**2 + (point_on_line.y-table['b'][i2])**2)**0.5)
    
    return total_offset

molecules = ('13CO_2014_merge', 'C18O_2014_merge', 'H2CO_303_202_bl', 'SiO_54_bl')

filenames = [molpath('APEX_{0}.fits'.format(molecule))
            for molecule in molecules]

molecules = molecules + ('H2CO_TemperatureFromRatio',
                         'H2CO_TemperatureFromRatio_smooth',
                         'H2CO_DendrogramTemperature',
                         'H2CO_DendrogramTemperature_smooth')
filenames.append(hpath('TemperatureCube_PiecewiseFromRatio.fits'))
filenames.append(hpath('TemperatureCube_smooth_PiecewiseFromRatio.fits'))
filenames.append(hpath('TemperatureCube_DendrogramObjects_Piecewise.fits'))
filenames.append(hpath('TemperatureCube_DendrogramObjects_smooth_Piecewise.fits'))


cmap = copy.copy(pl.cm.RdYlBu_r)
cmap.set_bad((1.0,)*3)
cmap.set_under((0.9,0.9,0.9,0.5))

for molecule,fn in zip(molecules,filenames):
    log.info(molecule)
    cube = spectral_cube.SpectralCube.read(fn)

    # respect_nan = False so that the background is zeros where there is data
    # and nan where there is not data
    # But not respecting nan results in data getting averaged with nan, so we
    # need to respect it and then manually flag (in a rather unreliable
    # fashion!)
    #pv = pvextractor.extract_pv_slice(cube, P, respect_nan=False)
    pv = pvextractor.extract_pv_slice(cube, P, respect_nan=True)
    bad_cols = np.isnan(np.nanmax(pv.data, axis=0))
    nandata = np.isnan(pv.data)
    pv.data[nandata & ~bad_cols] = 0

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

    for color, segment in zip(('red','green','blue','black','purple'),
                              ('abcde')):
        selection = table['segment'] == segment
        # Connect the previous to the next segment - force continuity
        if np.argmax(selection) > 0:
            selection[np.argmax(selection)-1] = True
        F.show_lines(np.array([[cdist[selection],
                                table["v'los"][selection]*1e3]]), zorder=1000,
                     color=color, linewidth=3, alpha=0.25)
        #F.show_markers(cdist[selection], table["v'los"][selection]*1e3, zorder=1000,
        #             color=color, marker='+')
    F.recenter(x=4.5/2., y=0., width=4.5, height=300000)
    #F.show_markers([offset_to_point(0.47,-0.01)], [30.404e3], color=['r'])
    #F.show_markers([offset_to_point(0.38,+0.04)], [39.195e3], color=['b'])
    F.show_markers([offset_to_point(0.47, -0.01)],[30.404e3], edgecolor='r', marker='x')
    F.show_markers([offset_to_point(0.38, 0.04)], [39.195e3], edgecolor='b', marker='x')
    F.show_markers([offset_to_point(0.253, 0.016)], [36.5e3], edgecolor='purple', marker='x')
    F.save(fpath('orbits/KDL2014_orbit_on_{0}.pdf'.format(molecule)))

    fig2 = pl.figure(2)
    pl.clf()
    img = cube.mean(axis=0).hdu
    img.data[np.isnan(img.data)] = 0
    F2 = aplpy.FITSFigure(img, convention='calabretta', figure=fig2)
    if 'Temperature' in fn:
        F2.show_colorscale(cmap=cmap, vmin=20, vmax=200)
        F2.add_colorbar()
    else:
        F2.show_grayscale()

    patches = P.to_patches(1, ec='gray', fc='none',
                           #transform=ax.transData,
                           clip_on=True, #clip_box=ax.bbox,
                           wcs=cube.wcs)
    for patch in patches:
        patch.set_linewidth(0.5)
        patch.set_alpha(0.2)
        patch.zorder = 50

    patches[0].set_edgecolor('green')
    patches[0].set_alpha(1)
    patches[0].set_linewidth(1)
    patches[0].zorder += 1

    patchcoll = matplotlib.collections.PatchCollection(patches, match_original=True)
    patchcoll.zorder=10

    c = F2._ax1.add_collection(patchcoll)

    F2._rectangle_counter += 1
    rectangle_set_name = 'rectangle_set_' + str(F2._rectangle_counter)

    for color, segment in zip(('red','green','blue','black','purple'),
                              ('abcde')):
        selection = table['segment'] == segment
        # Connect the previous to the next segment - force continuity
        if np.argmax(selection) > 0:
            selection[np.argmax(selection)-1] = True
        F2.show_lines(np.array([[table['l'][selection],
                                 table["b"][selection]]]), zorder=1000,
                     color=color, linewidth=3, alpha=0.5)
        #F2.show_markers(table['l'][selection], table["b"][selection], zorder=1000,
        #             color=color, marker='+', alpha=0.5)

    F2._layers[rectangle_set_name] = c
    F2.recenter(0, -0.03, width=1.8, height=0.3)
    F2.set_tick_labels_format('d.dd','d.dd')

    F2.show_markers([0.47], [-0.01], edgecolor='r', marker='x', zorder=1500)
    F2.show_markers([0.38], [0.04], edgecolor='b', marker='x', zorder=1500)
    F2.show_markers([0.253], [0.016], edgecolor='purple', marker='x', zorder=1500)

    F2.save(fpath('orbits/KDL2014_orbitpath_on_{0}.pdf'.format(molecule)))

    # Compute the temperature as a function of time in a ppv tube
    offset = np.linspace(0, cdist.max(), pv.shape[1])
    time = np.interp(offset, cdist, table['t'])
    vel = np.interp(time, table['t'], table["v'los"])
    y,x = np.indices(pv.data.shape)
    p,v = WCS(pv.header).wcs_pix2world(x,y, 0)

    vdiff = 15

    velsel = (v > (vel-vdiff)*1e3) & (v < (vel+vdiff)*1e3)

    pv.data[nandata] = np.nan
    pv.data[pv.data==0] = np.nan
    pv.data[~velsel] = np.nan

    mean_tem = np.nanmean(pv.data, axis=0)
    min_tem = np.nanmin(pv.data, axis=0)
    max_tem = np.nanmax(pv.data, axis=0)
    std_tem = np.nanstd(pv.data, axis=0)

    # errorbar version: ugly
    #eb = ax3.errorbar(time, mean_tem, yerr=[min_tem, max_tem],
    #                  linestyle='none', capsize=0, color='r', errorevery=20)
    #eb[-1][0].set_linestyle('--')

    #ax3.fill_between(time, mean_tem-std_tem, mean_tem+std_tem,
    #                 color='b', alpha=0.2)

    #ax3.plot(time, mean_tem, color='b', alpha=0.2)

    fig3 = pl.figure(3)
    fig3.clf()
    ax3 = fig3.gca()

    reftime = -2
    bricktime = 0.3
    ax3.plot(time-reftime, pv.data.T, 'k.', alpha=0.5, markersize=3)
    ax3.set_xlabel("Time since 1$^\\mathrm{st}$ pericenter passage [Myr]", size=24, labelpad=10)
    if 'Temperature' in fn:
        ax3.set_ylim(0,150)
        ax3.set_ylabel("Temperature [K]", size=24, labelpad=10)
        ytext = 135
    else:
        ax3.set_ylabel("$T_A^*$ [K]", size=24, labelpad=10)
        ytext = ax3.get_ylim()[1]*(14./15.)
    ax3.text(bricktime, ytext, "Brick", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='purple', weight='bold')
    ax3.text(bricktime+0.43, ytext, "Sgr B2", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='b', weight='bold')
    ax3.text(bricktime+3.58, ytext*135./140., "20 km s$^{-1}$", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='g', weight='bold')
    ax3.text(bricktime+3.66, ytext*135./140., "50 km s$^{-1}$", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='g', weight='bold')
    ax3.text(bricktime+3.28, ytext, "Sgr C", verticalalignment='center',
             horizontalalignment='center', rotation='vertical', color='m', weight='bold')
    pl.setp(ax3.get_xticklabels(), fontsize=20)
    pl.setp(ax3.get_yticklabels(), fontsize=20)
    ax3.set_xlim(-0.1,4.6)
    fig3.savefig(fpath('orbits/KDL2014_{0}_vs_time.pdf'.format(molecule)),
                 bbox_inches='tight')
