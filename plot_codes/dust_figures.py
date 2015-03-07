import numpy as np
import pyregion
from astropy import units as u
from paths import h2copath,figurepath,hpath,rpath,fpath,dpath
import copy
import os
import aplpy
import pylab as pl
import matplotlib
from astropy.io import fits
from masked_cubes import (cube303, cube303m, cube321, cube321m, cube303sm,
                          cube303msm, cube321sm, cube321msm)
from astropy import visualization
from astropy.visualization import mpl_normalize
from astropy import log
import wcsaxes
from aplpy.regions import ds9
from aplpy import regions as aplpyregions
from astropy.wcs import WCS
from agpy import asinh_norm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits import axisartist as AA

small_recen = dict(x=0.3, y=-0.03,width=1.05,height=0.27)
big_recen = dict(x=0.55, y=-0.075,width=2.3,height=0.40)

pl.ioff()

cm = matplotlib.cm.RdYlBu_r
cm.set_bad('#888888')
figsize=(20,10)

dustcolumn = fits.open('/Users/adam/work/gc/gcmosaic_column_conv36.fits')
dustcoldata = dustcolumn[0].data
dustcolwcs = WCS(dustcolumn[0].header)
dusttemperature = fits.open('/Users/adam/work/gc/gcmosaic_temp_conv36.fits')
dusttemdata = dusttemperature[0].data
dusttemwcs = WCS(dusttemperature[0].header)

gastemperature = fits.open(dpath('H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e4_masked.fits',))
gastemdata = gastemperature[0].data
gastemwcs = WCS(gastemperature[0].header)

for stuff in zip((dustcoldata, dusttemdata, dusttemdata, gastemdata),
                 (dustcolwcs,dusttemwcs,dusttemwcs, gastemwcs),
                 ('Column','Temperature', 'Temperature', 'Temperature',),
                 ('Column','Temperature', 'Temperature_FullRange', 'GasTemperature'),
                 (asinh_norm.AsinhNorm, matplotlib.colors.Normalize,
                  matplotlib.colors.Normalize, matplotlib.colors.Normalize),
                 ((0,100),(15,45),(10,200),(10,200))):

    im, wcs, label, fn, stretch, (vmin,vmax) = stuff

    fig1 = pl.figure(1, figsize=figsize)
    fig1.clf()
    ax = wcsaxes.WCSAxesSubplot(fig1, 1,1,1, wcs=wcs)
    fig1.add_axes(ax)
    norm = mpl_normalize.ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch())
    ims = ax.imshow(im, cmap=cm, vmin=vmin, vmax=vmax, norm=norm)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05,
                              axes_class=matplotlib.axes.Axes)

    cb = pl.colorbar(mappable=ims, cax=cax)
    cax.set_ylabel(label)
    cb.ax.yaxis.set_label_position('right')

    ax.set_xlabel("Galactic Longitude")
    ax.set_ylabel("Galactic Latitude")

    (xlo, ylo),(xhi,yhi) = wcs.wcs_world2pix([(big_recen['x']+big_recen['width']/2.,
                                               big_recen['y']-big_recen['height']/2.),
                                              (big_recen['x']-big_recen['width']/2.,
                                               big_recen['y']+big_recen['height']/2.)],
                                              0)
    ax.axis([xlo,xhi,ylo,yhi])
    fig1.savefig(fpath("big_maps/dust_{0}_color.png".format(fn)), bbox_inches='tight')

    (xlo, ylo),(xhi,yhi) = wcs.wcs_world2pix([(small_recen['x']+small_recen['width']/2.,
                                               small_recen['y']-small_recen['height']/2.),
                                              (small_recen['x']-small_recen['width']/2.,
                                               small_recen['y']+small_recen['height']/2.)],
                                              0)
    ax.axis([xlo,xhi,ylo,yhi])
    fig1.savefig(fpath("big_maps/dust_{0}_color_zoom.png".format(fn)), bbox_inches='tight')
