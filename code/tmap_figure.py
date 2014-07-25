import pylab as pl
import aplpy
import os
import copy
from paths import h2copath, figurepath

fig = pl.figure(4, figsize=(12,6))
fig.clf()
F = aplpy.FITSFigure(h2copath+'H2CO_321220_to_303202_smooth_bl_integ_temperature.fits',
                     convention='calabretta',
                     figure=fig)

cm = copy.copy(pl.cm.rainbow)
cm.set_bad('#888888')
F.show_colorscale(cmap=cm,vmin=15,vmax=200)
F.set_tick_labels_format('d.dd','d.dd')
F.recenter(0.3,-0.03,width=1.2,height=0.30)
peaksn = os.path.join(h2copath,'APEX_H2CO_303_202_smooth_mask_integ.fits')
F.show_contour(peaksn, levels=[4,7,11,20,38], colors=[(0,0,0,0.5)]*5, smooth=3,
               zorder=10, convention='calabretta')
F.add_colorbar()
F.colorbar.set_axis_label_text('T (K)')
F.save(os.path.join(figurepath, 'lores_tmap_withcontours.png'))
F.hide_layer('contour_set_1')
dusttemperature = '/Users/adam/work/gc/gcmosaic_temp_conv25.fits'
F.show_contour(dusttemperature,
               levels=[20,25],
               colors=[(0,0,x,0.5) for x in [0.8,0.6,0.4,0.2]], zorder=10)
F.save(os.path.join(figurepath,'lores_tmap_withtdustcontours.png'))
#F.show_contour('h2co218222_all.fits', levels=[1,7,11,20,38], colors=['g']*5, smooth=1, zorder=5)
#F.show_contour(datapath+'APEX_H2CO_merge_high_smooth_noise.fits', levels=[0.05,0.1], colors=['#0000FF']*2, zorder=3, convention='calabretta')
#F.show_contour(datapath+'APEX_H2CO_merge_high_nhits.fits', levels=[9], colors=['#0000FF']*2, zorder=3, convention='calabretta',smooth=3)
#F.show_regions('2014_expansion_targets_simpler.reg')
#F.save('CMZ_H2CO_observed_planned.png')
#F.show_rgb(background, wcs=wcs)
#F.save('CMZ_H2CO_observed_planned_colorful.png')


fig = pl.figure(5, figsize=(12,6))
fig.clf()
F2 = aplpy.FITSFigure(dusttemperature, convention='calabretta', figure=fig)
F2.show_colorscale(cmap=pl.cm.hot, vmin=10, vmax=40)
F2.add_colorbar()
F2.show_contour(h2copath+'H2CO_321220_to_303202_smooth_bl_integ_temperature.fits',
                convention='calabretta',
                levels=[30,75,100,150],
                cmap=pl.cm.BuGn)
F2.recenter(0.3,-0.03,width=1.2,height=0.30)

F2.save(os.path.join(figurepath,'H2COtemperatureOnDust.png'))
