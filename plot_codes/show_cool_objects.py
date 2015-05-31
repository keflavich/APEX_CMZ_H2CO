import pylab as pl
import numpy as np
import aplpy
import os
import copy
from astropy import log
from paths import h2copath, figurepath
from astropy.table import Table
import paths
import matplotlib
from scipy import stats as ss
from astropy.io import fits
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

pl.ioff()
# Close these figures so we can remake them in the appropriate size
for fignum in (4,5,6,7):
    pl.close(fignum)

cmap = pl.cm.RdYlBu_r
figsize = (20,10)

small_recen = dict(x=0.3, y=-0.03,width=1.05,height=0.27)
big_recen = dict(x=0.55, y=-0.075,width=2.3,height=0.40)

sgrb2x = [000.6773, 0.6578, 0.6672]
sgrb2y = [-00.0290, -00.0418, -00.0364]

vmin=10
vmax = 200

ftemplate = 'H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e4_masked.fits'

fig = pl.figure(4, figsize=figsize)
fig.clf()
smooth =''
F = aplpy.FITSFigure(h2copath+ftemplate.format(smooth),
                     convention='calabretta',
                     figure=fig)

cm = copy.copy(cmap)
cm.set_bad((0.5,)*3)
F.show_colorscale(cmap=cm,vmin=vmin,vmax=vmax)
F.set_tick_labels_format('d.dd','d.dd')

tbl = Table.read(paths.tpath('PPV_H2CO_Temperature.ipac'), format='ascii.ipac')
cooler_than_dust = ((tbl['higaldusttem'] > tbl['temperature_chi2']) &
                    (tbl['IsNotH2CO'] == 'False') & (tbl['is_leaf']=='True') &
                    (tbl['gausscorrfactor'] < 3))
cooler_than_turb = ((tbl['DespoticTem'] < tbl['tmin1sig_chi2']) &
                    (tbl['IsNotH2CO'] == 'False') & (tbl['is_leaf']=='True') &
                    (tbl['gausscorrfactor'] < 3))
F.show_markers(tbl[cooler_than_dust]['x_cen']-360*(tbl[cooler_than_dust]['x_cen']>180),
               tbl[cooler_than_dust]['y_cen'], marker='x', color='k',
               facecolor='k', edgecolor='k')
F.show_markers(tbl[cooler_than_turb]['x_cen']-360*(tbl[cooler_than_turb]['x_cen']>180),
               tbl[cooler_than_turb]['y_cen'], marker='+', color='k',
               facecolor='k', edgecolor='k')


# PV
recen = dict(x=0.57, y=20e3, width=2.26, height=250e3)
vmin=10
vmax = 200
cmap = pl.cm.RdYlBu_r
figsize = (20,10)
ftemplate = 'pv_H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e4.fits'
hdu = fits.open(paths.hpath(ftemplate.format(smooth)))[0]
hdu.header['CTYPE1'] = 'GLON'

fig = pl.figure(5, figsize=figsize)
fig.clf()

F = aplpy.FITSFigure(hdu,
                     figure=fig)

cm = copy.copy(cmap)
cm.set_bad((0.5,)*3)
F.show_colorscale(cmap=cm,vmin=vmin,vmax=vmax,
                  aspect=1 if smooth=='' else 2)
F.recenter(**recen)
F.tick_labels.set_xformat('d.dd')
F.add_colorbar()
F.colorbar.set_ticks(np.arange(20,240,40))
F.colorbar.set_axis_label_font(size=20)
F.colorbar.set_axis_label_text('Temperature (K)')

F.show_markers(tbl[cooler_than_dust]['x_cen']-360*(tbl[cooler_than_dust]['x_cen']>180),
               tbl[cooler_than_dust]['v_cen'], marker='x', color='k',
               facecolor='k', edgecolor='k')
F.show_markers(tbl[cooler_than_turb]['x_cen']-360*(tbl[cooler_than_turb]['x_cen']>180),
               tbl[cooler_than_turb]['v_cen'], marker='+', color='k',
               facecolor='k', edgecolor='k')

F.refresh()
F._ax1.set_ylabel("$V_{LSR} (\mathrm{km\ s}^{-1})$")
#F._ax1.set_yticklabels([(str(int(x.get_text())/1000)) for x in F._ax1.get_yticklabels()])
F.refresh()
F.axis_labels.set_xtext(r'Galactic Longitude')

pl.draw()
pl.show()
