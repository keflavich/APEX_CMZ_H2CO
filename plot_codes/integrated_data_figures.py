import numpy as np
import pyregion
from astropy import units as u
from paths import h2copath,figurepath,hpath,rpath,fpath
import copy
import os
import aplpy
import pylab as pl
import matplotlib
from astropy.io import fits
from masked_cubes import (cube303, cube303m, cube321, cube321m, cube303sm,
                          cube303msm, cube321sm, cube321msm)
from astropy import log
import wcsaxes
from aplpy.regions import ds9
from aplpy import regions as aplpyregions
from astropy.wcs import WCS
from agpy import asinh_norm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits import axisartist as AA


pl.ioff()

cm = matplotlib.cm.RdYlBu_r
cm.set_bad('#888888')
figsize=(20,10)

method_label = {'mean': "Mean Brightness $T_A^*$ (K)",
                'moment0': "Integrated Brightness\n"
                   r"$T_A^* dv $ (K km s$^{-1}$)",
                'max': 'Peak Brightness $T_A^*$ (K)'}

dustcolumn = fits.open('/Users/adam/work/gc/gcmosaic_column_conv36.fits')
dustcoldata = dustcolumn[0].data
dustcolwcs = WCS(dustcolumn[0].header)

# Direct from cubes
for method in ('max','mean','moment0'):
    for cube,name,masked,smooth in zip((cube303, cube303m, cube303sm, cube303msm,
                                        cube321, cube321m,  cube321sm, cube321msm),
                                       ['303_202']*4 + ['321_220']*4,
                                       ['','_masked']*4,
                                       ['','','_smooth','_smooth']*2):

        outname = "APEX_H2CO_{0}{1}{2}_{3}.fits".format(name,masked,smooth,method)
        log.info(outname)
        if not os.path.exists(hpath(outname)):
            im = getattr(cube,method)(axis=0).hdu
            im.writeto(hpath(outname), clobber=True)
            log.info("{0} written.".format(outname))
        else:
            im = fits.open(hpath(outname))[0]
            log.info("{0} read.".format(outname))

        wcs = WCS(im.header)

        fig1 = pl.figure(1, figsize=figsize)
        fig1.clf()
        ax = wcsaxes.WCSAxesSubplot(fig1, 1,1,1, wcs=wcs)
        fig1.add_axes(ax)

        cm = pl.cm.gray_r
        #im.data -= (np.nanmin(im.data) - 0.01)
        vmin = np.percentile(im.data[np.isfinite(im.data)], 00.05)
        vmin = 0
        vmax = np.percentile(im.data[np.isfinite(im.data)], 99.95)
        ims = ax.imshow(im.data, cmap=cm, vmin=vmin, vmax=vmax)
                        #norm=matplotlib.colors.LogNorm())
        limits = ax.axis()


        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        # http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05,
                                  axes_class=matplotlib.axes.Axes)

        cb = pl.colorbar(mappable=ims, cax=cax)
        if len(cax.yaxis.get_ticklabels()) < 2:
            ticks = np.logspace(np.log10(vmin), np.log10(vmax), 5)
            cax.yaxis.set_ticks(ticks)
            cax.yaxis.set_ticklabels(["{0:0.2g}".format(x) for x in ticks])
        cax.set_ylabel(method_label[method])
        cb.ax.yaxis.set_label_position('right')

        ax.set_xlabel("Galactic Longitude")
        ax.set_ylabel("Galactic Latitude")

        tr = ax.get_transform(dustcolwcs)
        con = ax.contour(dustcoldata, levels=[5], colors=[(0,0,0,0.5)],
                         zorder=15, alpha=0.5, linewidths=[0.5], transform=tr)
        ax.axis(limits)

        pl.savefig(fpath('integrated/{0}'.format(outname.replace(".fits",".pdf"))),
                   bbox_inches='tight')
        for c in con.collections: c.set_visible(False)
        
        labels = pyregion.open(rpath('ridge_names.reg'))
        PC, TC = ds9(labels, im.header, text_offset=0)
        #PC.add_to_axes(ax)
        TC = aplpyregions.ArtistCollection([x for x in TC.artistlist if isinstance(x, matplotlib.text.Annotation)])
        TC.add_to_axes(ax)

        pl.savefig(fpath('integrated/{0}'.format(outname.replace(".fits","_labeled.pdf"))),
                   bbox_inches='tight')


        for c in con.collections: c.set_visible(True)
        ax.axis(limits)
        pl.savefig(fpath('integrated/{0}'.format(outname.replace(".fits",
                                                                 "_labeled_dustcontours.pdf"))),
                   bbox_inches='tight')


# From make_ratiotem_cubesims results
for smooth in ("","_bl", "_smooth", "_smooth_bl"):

    outname = "APEX_H2CO_303_202{0}_mask_integ.fits".format(smooth)
    im = fits.open(hpath(outname))[0]
    log.info("{0} read.".format(outname))

    wcs = WCS(im.header)

    fig1 = pl.figure(1, figsize=figsize)
    fig1.clf()
    ax = wcsaxes.WCSAxesSubplot(fig1, 1,1,1, wcs=wcs)
    fig1.add_axes(ax)

    cm = pl.cm.gray_r
    #im.data -= np.nanmin(im.data) - 0.00001
    vmin = np.percentile(im.data[np.isfinite(im.data)], 00.05)
    vmin = 0
    vmax = np.percentile(im.data[np.isfinite(im.data)], 99.95)
    ims = ax.imshow(im.data, cmap=cm, vmin=vmin, vmax=vmax)
                    #norm=matplotlib.colors.LogNorm())


    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    # http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="5%", pad=0.05, axes_class=matplotlib.axes.Axes)
    fig1.add_axes(cax)
    cax.yaxis.tick_right()

    cb = pl.colorbar(mappable=ims, cax=cax)
    #cb.set_ticks(np.linspace(vmin, vmax, 5))
    cax.set_ylabel(method_label['moment0'])
    cax.yaxis.set_label_position('right')
    for tl in cax.get_xticklabels():
        tl.set_visible(False)
    cax.yaxis.tick_right()

    ax.set_xlabel("Galactic Longitude")
    ax.set_ylabel("Galactic Latitude")

    pl.savefig(fpath('integrated/{0}'.format(outname.replace("fits","pdf"))),
               bbox_inches='tight')

    tr = ax.get_transform(dustcolwcs)
    con = ax.contour(dustcoldata, levels=[5], colors=[(0,0,0,0.5)],
                     zorder=15, alpha=0.5, linewidths=[0.5], transform=tr)

    pl.savefig(fpath('integrated/{0}'.format(outname.replace(".fits",
                                                             "_dustcontours.pdf"))),
               bbox_inches='tight')
