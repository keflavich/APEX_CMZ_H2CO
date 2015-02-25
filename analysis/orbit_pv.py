import os
import numpy as np
import pvextractor
from pvextractor.geometry.poly_slices import extract_poly_slice
import spectral_cube
from spectral_cube import SpectralCube, BooleanArrayMask, Projection
import aplpy
import pylab as pl
import matplotlib
import copy
from paths import mpath,apath,fpath,molpath,hpath
from astropy import units as u
from astropy import coordinates
from astropy.io import ascii, fits
from astropy import log
from astropy.wcs import WCS
from astropy import wcs
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

reftime = -2
bricktime = 0.3
time = table['t']

# how much time per pixel?
dtdx = (table['t'].max() - table['t'].min()) / cdist.max()

figsize=(20,10)

vmin,vmax=10,200

def offset_to_point(glon, glat):
    """
    Determine the offset along the orbit to the nearest point on an orbit to
    the specified point
    """
    import shapely.geometry as geom

    line = geom.LineString(table['l','b'])
    point = geom.Point(glon, glat)
    return line.project(point)

molecules = ('13CO_2014_merge', 'C18O_2014_merge', 'H2CO_303_202_bl', 'SiO_54_bl')

filenames = [molpath('APEX_{0}.fits'.format(molecule))
            for molecule in molecules]

molecules = molecules + ('H2CO_TemperatureFromRatio',
                         'H2CO_TemperatureFromRatio_smooth',
                         'H2CO_Ratio',
                         'H2CO_Ratio_smooth',
                         #'H2CO_DendrogramTemperature',
                         #'H2CO_DendrogramTemperature_smooth',
                         #'H2CO_DendrogramTemperature_Leaves',
                         #'H2CO_DendrogramTemperature_Leaves_smooth'
                        )
filenames.append(hpath('TemperatureCube_PiecewiseFromRatio.fits'))
filenames.append(hpath('TemperatureCube_smooth_PiecewiseFromRatio.fits'))
filenames.append(hpath('H2CO_321220_to_303202_cube_bl.fits'))
filenames.append(hpath('H2CO_321220_to_303202_cube_smooth_bl.fits'))
#filenames.append(hpath('TemperatureCube_DendrogramObjects.fits'))
#filenames.append(hpath('TemperatureCube_DendrogramObjects_smooth.fits'))
#filenames.append(hpath('TemperatureCube_DendrogramObjects_leaves.fits'))
#filenames.append(hpath('TemperatureCube_DendrogramObjects_smooth_leaves.fits'))


cmap = copy.copy(pl.cm.RdYlBu_r)
cmap.set_under((0.9,0.9,0.9,0.5))

for weight in ("_weighted",""):
    for molecule,fn in zip(molecules[-4:],filenames[-4:]):
        log.info(molecule)
        cube = spectral_cube.SpectralCube.read(fn)

        if 'smooth' in fn:
            cmap.set_bad((1.0,)*3)
        else:
            cmap.set_bad((0.9,0.9,0.9,0.5))


        if weight:
            wcube = (spectral_cube.SpectralCube.read(hpath('APEX_H2CO_303_202_smooth_bl.fits'))
                     if 'smooth' in fn
                     else spectral_cube.SpectralCube.read(hpath('APEX_H2CO_303_202_bl.fits')))
            mcube = (SpectralCube.read(hpath('APEX_H2CO_303_202_smooth_bl_mask.fits'))
                     if 'smooth' in fn
                     else SpectralCube.read(hpath('APEX_H2CO_303_202_bl_mask.fits')))
            if cube.shape != wcube.shape:
                log.info("Not weighting {0}".format(fn))
                continue
            weighted = copy.copy(cube)
            weighted._data = wcube._data * cube._data
            mask = (wcube.mask & cube.mask &
                    BooleanArrayMask(weighted.filled_data[...] != 0, weighted.wcs) &
                    BooleanArrayMask(wcube.filled_data[...] != 0, wcube.wcs) &
                    BooleanArrayMask(mcube._data==1, mcube.wcs)
                   )
            weighted = weighted.with_mask(mask)
            wcube = wcube.with_mask(mask)

        pvfilename = paths.dpath('orbits/KDL2014_orbit_on_{0}{1}.fits'.format(molecule, weight))
        if os.path.exists(pvfilename):
            pv = fits.open(pvfilename)[0]
        else:
            if weight:

                #pv1a = pvextractor.extract_pv_slice(cube, P, respect_nan=True)
                pv1,apv1 = extract_poly_slice(weighted.filled_data[...].value,
                                              P.sample_polygons(spacing=1.0,
                                                                wcs=cube.wcs),
                                              return_area=True)
                print()
                pv2,apv2 = extract_poly_slice(wcube.filled_data[...].value,
                                         P.sample_polygons(spacing=1.0, wcs=cube.wcs),
                                         return_area=True)

                assert np.count_nonzero(pv1) == np.count_nonzero(pv2)
                assert np.all((pv1 != 0) == (pv2 != 0))
                assert np.all(apv1==apv2)

                # Some of this header manipulation seems unnecessary but isn't.
                header = copy.copy(cube.header)
                for kw in ('CRPIX','CRVAL','CDELT','CUNIT','CTYPE',):
                    header[kw+"2"] = header[kw+"3"]
                    del header[kw+"3"]
                header['CTYPE1'] = 'OFFSET'
                header = pvextractor.utils.wcs_slicing.slice_wcs(WCS(header),
                                                                 spatial_scale=7.2*u.arcsec).to_header()
                pv = fits.PrimaryHDU(data=pv1/pv2,
                                     header=header)
                pv2hdu = fits.PrimaryHDU(data=pv2, header=header)
            else:
                # respect_nan = False so that the background is zeros where there is data
                # and nan where there is not data
                # But not respecting nan results in data getting averaged with nan, so we
                # need to respect it and then manually flag (in a rather unreliable
                # fashion!)
                #pv = pvextractor.extract_pv_slice(cube, P, respect_nan=False)
                pv = pvextractor.extract_pv_slice(cube, P, respect_nan=True)
            if not os.path.isdir(os.path.dirname(pvfilename)):
                os.mkdir(os.path.dirname(pvfilename))
            pv.writeto(pvfilename)
        bad_cols = np.isnan(np.nanmax(pv.data, axis=0))
        nandata = np.isnan(pv.data)
        pv.data[nandata & ~bad_cols] = 0

        fig1 = pl.figure(1, figsize=figsize)
        fig1.clf()
        ax = fig1.gca()
        mywcs = WCS(pv.header)
        xext, = mywcs.sub([1]).wcs_pix2world((0,pv.shape[1]), 0)
        yext, = mywcs.sub([2]).wcs_pix2world((0,pv.shape[0]), 0)
        yext /= 1e3
        dy = yext[1]-yext[0]
        dx = xext[1]-xext[0]
        #F = aplpy.FITSFigure(pv, figure=fig1)
        actual_aspect = pv.shape[0]/float(pv.shape[1])
        if 'Temperature' in fn:
            #F.show_colorscale(cmap=cmap, aspect=0.5/actual_aspect, vmin=vmin, vmax=vmax)
            im = ax.imshow(pv.data, extent=[xext[0], xext[1], yext[0], yext[1]],
                           aspect=0.5*dx/dy, vmin=vmin, vmax=vmax, cmap=cmap)
            #divider = make_axes_locatable(ax)
            #cax = divider.append_axes("right", size="2%", pad=0.05)
            cb = fig1.colorbar(im)
            cb.set_label("Temperature (K)")
            ## This is where it fails...
            ##F.add_colorbar()
            #if 'smooth' in fn:
            #    F.colorbar.set_pad(-0.94/actual_aspect)
            #else:
            #    F.colorbar.set_pad(-1.43/actual_aspect)
            #F.colorbar.set_axis_label_text("Temperature (K)")
            #pl.colorbar(F._ax1.images[0], cax=cax)
        elif 'Ratio' in fn:
            im = ax.imshow(pv.data, extent=[xext[0], xext[1], yext[0], yext[1]],
                           aspect=0.5*dx/dy, vmin=0, vmax=0.5, cmap=cmap)
            cb = fig1.colorbar(im)
            cb.set_label("Ratio $R_1$")
        else:
            #F.show_grayscale(aspect=0.5/actual_aspect)
            im = ax.imshow(pv.data, extent=[xext[0], xext[1], yext[0], yext[1]],
                           aspect=0.5*dx/dy, cmap=cmap)

        ax2 = ax.twiny()
        ax2.xaxis.set_label_position('top')
        ax2.xaxis.set_ticklabels(["{0:0.2f}".format(x) for x in
                                  np.interp(ax2.xaxis.get_ticklocs(),
                                            np.linspace(0,1,time.size),
                                            time-reftime)])
        ax2.set_xlabel("Time since 1$^\\mathrm{st}$ pericenter passage [Myr]",
                       size=24, labelpad=10)
        ax.set_xlabel("Offset (degrees)")
        ax.set_ylabel("$V_{LSR}$ $(\mathrm{km\ s}^{-1})$")

        for color, segment in zip(('red','green','blue','black','purple'),
                                  ('abcde')):
            selection = table['segment'] == segment
            # Connect the previous to the next segment - force continuity
            if np.argmax(selection) > 0:
                selection[np.argmax(selection)-1] = True
            ax.plot(np.array(cdist[selection]),
                    table["v'los"][selection], zorder=1000,
                    color=color, linewidth=3, alpha=0.25)
            #F.show_lines(np.array([[cdist[selection],
            #                        table["v'los"][selection]*1e3]]), zorder=1000,
            #             color=color, linewidth=3, alpha=0.25)
            #F.show_markers(cdist[selection], table["v'los"][selection]*1e3, zorder=1000,
            #             color=color, marker='+')
        #F.recenter(x=4.5/2., y=20., width=4.5, height=240000)
        #F.show_markers([offset_to_point(0.47,-0.01)], [30.404e3], color=['r'])
        #F.show_markers([offset_to_point(0.38,+0.04)], [39.195e3], color=['b'])
        #F.show_markers([offset_to_point(0.47, -0.01)],[30.404e3], edgecolor='r', marker='x')
        #F.show_markers([offset_to_point(0.38, 0.04)], [39.195e3], edgecolor='b', marker='x')
        #F.show_markers([offset_to_point(0.253, 0.016)], [36.5e3], edgecolor='purple', marker='x')
        ax.plot(offset_to_point(0.47, -0.01),30.404, color='r', marker='x')
        ax.plot(offset_to_point(0.38, 0.04), 39.195, color='b', marker='x')
        ax.plot(offset_to_point(0.253, 0.016), 36.5, color='purple', marker='x')

        #F.refresh()
        #F._ax1.set_ylabel("$V_{LSR} (\mathrm{km\ s}^{-1})$")
        #F._ax1.set_yticklabels([(str(int(x.get_text())/1000)) for x in F._ax1.get_yticklabels()])
        #F.refresh()
        ax.axis([xext[0], xext[1], yext[0], yext[1]])

        #F.save(fpath('orbits/KDL2014_orbit_on_{0}{1}.pdf'.format(molecule, weight)))
        fig1.savefig(fpath('orbits/KDL2014_orbit_on_{0}{1}.pdf'.format(molecule, weight)),
                     bbox_inches='tight')

        # TODO!
        #if weight:
        #    color = (0.5,)*3 # should be same as background #888
        #    ax.show_contour(pv2hdu,
        #                   levels=[-1,0]+np.logspace(0.20,2).tolist(),
        #                   colors=([(0.5,0.5,0.5,1)]*2 + 
        #                           [color + (alpha,) for alpha in
        #                            np.exp(-(np.logspace(0.20,2)-1.7)**2/(2.5**2*2.))]),
        #                   filled=True,
        #                   smooth=3,
        #                   zorder=10, convention='calabretta')
        #    F.save(fpath('orbits/KDL2014_orbit_on_{0}{1}_masked.pdf'.format(molecule, weight)))


        fig2 = pl.figure(2, figsize=figsize)
        pl.clf()
        if weight:
            im1 = weighted.sum(axis=0)
            im2 = wcube.sum(axis=0)
            im = im1.value / im2.value
            img = im1.hdu
            img.data = im

        else:
            img = cube.mean(axis=0).hdu
        img.data[np.isnan(img.data)] = 0
        F2 = aplpy.FITSFigure(img, convention='calabretta', figure=fig2)
        if 'Temperature' in fn:
            F2.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax)
            F2.add_colorbar()
            F2.colorbar.set_ticks(np.arange(20,240,40))
            F2.colorbar.set_axis_label_font(size=20)
            F2.colorbar.set_axis_label_text('Temperature (K)')
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

        F2.save(fpath('orbits/KDL2014_orbitpath_on_{0}{1}.pdf'.format(molecule, weight)))

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

        fig3 = pl.figure(3, figsize=figsize)
        fig3.clf()
        ax3 = fig3.gca()

        ax3.plot(time-reftime, pv.data.T, 'k.', alpha=0.5, markersize=3)
        ax3.set_xlabel("Time since 1$^\\mathrm{st}$ pericenter passage [Myr]", size=24, labelpad=10)
        if 'Temperature' in fn:
            ax3.set_ylim(0,vmax)
            ax3.set_ylabel("Temperature [K]", size=24, labelpad=10)
            ytext = 180
        elif 'Ratio' in fn:
            ax3.set_ylim(0,0.5)
            ax3.set_ylabel("Ratio $R_1$", size=24, labelpad=10)
            ytext = 0.5*14./15.
        else:
            ax3.set_ylabel("$T_A^*$ [K]", size=24, labelpad=10)
            ytext = ax3.get_ylim()[1]*(14./15.)

        ax3.text(bricktime, ytext, "Brick", verticalalignment='center',
                 horizontalalignment='center', rotation='vertical', color='k', weight='bold')
        ax3.text(bricktime+0.43, ytext, "Sgr B2", verticalalignment='center',
                 horizontalalignment='center', rotation='vertical', color='k', weight='bold')
        ax3.text(bricktime+3.58, ytext*135./140., "20 km s$^{-1}$", verticalalignment='center',
                 horizontalalignment='center', rotation='vertical', color='k', weight='bold')
        ax3.text(bricktime+3.66, ytext*135./140., "50 km s$^{-1}$", verticalalignment='center',
                 horizontalalignment='center', rotation='vertical', color='k', weight='bold')
        ax3.text(bricktime+3.28, ytext, "Sgr C", verticalalignment='center',
                 horizontalalignment='center', rotation='vertical', color='k', weight='bold')
        pl.setp(ax3.get_xticklabels(), fontsize=20)
        pl.setp(ax3.get_yticklabels(), fontsize=20)
        ax3.set_xlim(-0.1,4.6)
        fig3.savefig(fpath('orbits/KDL2014_{0}_vs_time_{1}.pdf'.format(molecule, weight)),
                     bbox_inches='tight')
