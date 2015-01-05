import pylab as pl
import numpy as np
import aplpy
import os
import copy
from astropy import log
from paths import h2copath, figurepath

# Close these figures so we can remake them in the appropriate size
for fignum in (4,5):
    pl.close(fignum)

for ftemplate,outtype in zip(('H2CO_321220_to_303202{0}_bl_integ_temperature.fits',
                              'TemperatureCube_DendrogramObjects{0}_Piecewise_integ.fits'),
                             ('','dendro')):

    for smooth in ("","_smooth",):#"_vsmooth"):
        fig = pl.figure(4, figsize=(14,7))
        fig.clf()
        F = aplpy.FITSFigure(h2copath+ftemplate.format(smooth),
                             convention='calabretta',
                             figure=fig)

        cm = copy.copy(pl.cm.rainbow)
        cm.set_bad((0.5,)*3)
        F.show_colorscale(cmap=cm,vmin=15,vmax=200)
        F.set_tick_labels_format('d.dd','d.dd')
        F.recenter(0.3,-0.03,width=1.2,height=0.30)
        peaksn = os.path.join(h2copath,'APEX_H2CO_303_202{0}_bl_mask_integ.fits'.format(smooth))
        #F.show_contour(peaksn, levels=[4,7,11,20,38], colors=[(0.25,0.25,0.25,0.5)]*5, #smooth=3,
        #               linewidths=[1.0]*5,
        #               zorder=10, convention='calabretta')
        #color = (0.25,)*3
        #F.show_contour(peaksn, levels=[4,7,11,20,38], colors=[color + (alpha,) for alpha in (0.9,0.6,0.3,0.1,0.0)], #smooth=3,
        #               filled=True,
        #               #linewidths=[1.0]*5,
        #               zorder=10, convention='calabretta')
        color = (0.5,)*3 # should be same as background #888
        F.show_contour(peaksn, levels=[0]+np.logspace(0.20,2).tolist(),
                       colors=[(0.5,0.5,0.5,1)] + [color + (alpha,) for alpha in np.exp(-(np.logspace(0.20,2)-1.7)**2/(2.5**2*2.))], #smooth=3,
                       filled=True,
                       #linewidths=[1.0]*5,
                       zorder=10, convention='calabretta')
        F.add_colorbar()
        F.colorbar.set_axis_label_text('T (K)')
        F.save(os.path.join(figurepath, "big_maps", 'lores{0}{1}_tmap_withcontours.pdf'.format(smooth, outtype)))
        F.recenter(0.55,-0.075,width=2.3,height=0.40)
        F.save(os.path.join(figurepath, "big_maps", 'big_lores{0}{1}_tmap_withcontours.pdf'.format(smooth, outtype)))
        log.info(os.path.join(figurepath, "big_maps", 'big_lores{0}{1}_tmap_withcontours.pdf'.format(smooth, outtype)))

        Fsn = aplpy.FITSFigure(peaksn, convention='calabretta')
        Fsn.show_grayscale(vmin=0, vmax=10, stretch='linear', invert=True)
        Fsn.add_colorbar()
        Fsn.colorbar.set_axis_label_text('Peak S/N')
        Fsn.set_tick_labels_format('d.dd','d.dd')
        Fsn.recenter(0.55,-0.075,width=2.3,height=0.40)
        Fsn.save(os.path.join(figurepath, "big_maps", 'big_lores{0}{1}_peaksn.pdf'.format(smooth, outtype)))


        F.hide_layer('contour_set_1')
        dusttemperature = '/Users/adam/work/gc/gcmosaic_temp_conv25.fits'
        F.show_contour(dusttemperature,
                       levels=[20,25],
                       colors=[(0,0,x,0.5) for x in [0.9,0.7,0.6,0.2]], zorder=10)
        F.recenter(0.3,-0.03,width=1.2,height=0.30)
        F.save(os.path.join(figurepath, "big_maps",'lores{0}{1}_tmap_withtdustcontours.pdf'.format(smooth, outtype)))
        F.recenter(0.55,-0.075,width=2.3,height=0.40)
        F.save(os.path.join(figurepath, "big_maps",'big_lores{0}{1}_tmap_withtdustcontours.pdf'.format(smooth, outtype)))
        log.info(os.path.join(figurepath, "big_maps",'big_lores{0}{1}_tmap_withtdustcontours.pdf'.format(smooth, outtype)))
    #F.show_contour('h2co218222_all.fits', levels=[1,7,11,20,38], colors=['g']*5, smooth=1, zorder=5)
    #F.show_contour(datapath+'APEX_H2CO_merge_high_smooth_noise.fits', levels=[0.05,0.1], colors=['#0000FF']*2, zorder=3, convention='calabretta')
    #F.show_contour(datapath+'APEX_H2CO_merge_high_nhits.fits', levels=[9], colors=['#0000FF']*2, zorder=3, convention='calabretta',smooth=3)
    #F.show_regions('2014_expansion_targets_simpler.reg')
    #F.save('CMZ_H2CO_observed_planned.pdf')
    #F.show_rgb(background, wcs=wcs)
    #F.save('CMZ_H2CO_observed_planned_colorful.pdf')


fig = pl.figure(5, figsize=(14,7))
fig.clf()
F2 = aplpy.FITSFigure(dusttemperature, convention='calabretta', figure=fig)
F2.show_colorscale(cmap=pl.cm.hot, vmin=10, vmax=40)
F2.add_colorbar()
F2.show_contour(h2copath+'H2CO_321220_to_303202_smooth_bl_integ_temperature.fits',
                convention='calabretta',
                levels=[30,75,100,150],
                cmap=pl.cm.BuGn)
F2.recenter(0.3,-0.03,width=1.2,height=0.30)

F2.save(os.path.join(figurepath, "big_maps",'H2COtemperatureOnDust.pdf'))

F2.recenter(0.55,-0.075,width=2.3,height=0.40)
F2.save(os.path.join(figurepath, "big_maps",'big_H2COtemperatureOnDust.pdf'))

fig = pl.figure(6, figsize=(14,7))
fig.clf()
F = aplpy.FITSFigure('/Users/adam/work/gc/Tkin-GC.fits.gz',
                     convention='calabretta',
                     figure=fig)

cm = copy.copy(pl.cm.rainbow)
cm.set_bad((0.5,)*3)
F.show_colorscale(cmap=cm,vmin=15,vmax=200)
F.set_tick_labels_format('d.dd','d.dd')
F.recenter(0.3,-0.03,width=1.2,height=0.30)
F.add_colorbar()
F.colorbar.set_axis_label_text('T (K)')
F.save(os.path.join(figurepath, "big_maps", 'ott2014_nh3_tmap_15to200.pdf'))
F.show_colorscale(cmap=cm,vmin=15,vmax=80)
F.save(os.path.join(figurepath, "big_maps", 'ott2014_nh3_tmap_15to80.pdf'))
