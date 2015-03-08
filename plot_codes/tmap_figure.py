import pylab as pl
import numpy as np
import aplpy
import os
import copy
from astropy import log
from paths import h2copath, figurepath
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

dustcolumn = '/Users/adam/work/gc/gcmosaic_column_conv36.fits'

toloop = zip((
              'H2CO_321220_to_303202{0}_bl_integ_temperature_dens3e4.fits',
              'H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens3e4.fits',
              'H2CO_321220_to_303202{0}_bl_integ_temperature_dens1e4.fits',
              'H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e4.fits',
              'H2CO_321220_to_303202{0}_bl_integ_temperature_dens1e5.fits',
              'H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e5.fits',
              'H2CO_321220_to_303202{0}_bl_integ_temperature_dens1e4_masked.fits',
              'H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e4_masked.fits',
              'H2CO_321220_to_303202{0}_bl_integ_temperature_dens3e4_masked.fits',
              'H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens3e4_masked.fits',
              'H2CO_321220_to_303202{0}_bl_integ_temperature_dens1e5_masked.fits',
              'H2CO_321220_to_303202{0}_bl_integ_weighted_temperature_dens1e5_masked.fits',
              'TemperatureCube_DendrogramObjects{0}_leaves_integ.fits',
              'TemperatureCube_DendrogramObjects{0}_leaves_integ_weighted.fits',
              'TemperatureCube_DendrogramObjects{0}_integ.fits',
              'TemperatureCube_DendrogramObjects{0}_integ_weighted.fits'),
             ('dens3e4',       'dens3e4_weighted',
              'dens1e4',       'dens1e4_weighted',
              'dens1e5',       'dens1e5_weighted',
              'dens1e4_masked','dens1e4_weighted_masked',
              'dens3e4_masked','dens3e4_weighted_masked',
              'dens1e5_masked','dens1e5_weighted_masked',
              'dendro_leaf','dendro_leaf_weighted',
              'dendro','dendro_weighted'))


for vmax,vmax_str in zip((100,200),("_vmax100","")):
    for ftemplate,outtype in toloop:

        for smooth in ("","_smooth",):#"_vsmooth"):
            log.info(ftemplate.format(smooth)+"   "+outtype)
            fig = pl.figure(4, figsize=figsize)
            fig.clf()
            F = aplpy.FITSFigure(h2copath+ftemplate.format(smooth),
                                 convention='calabretta',
                                 figure=fig)

            cm = copy.copy(cmap)
            cm.set_bad((0.5,)*3)
            F.show_colorscale(cmap=cm,vmin=vmin,vmax=vmax)
            F.set_tick_labels_format('d.dd','d.dd')
            F.recenter(**small_recen)
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
            F.show_contour(peaksn, levels=[-1,0]+np.logspace(0.20,2).tolist(),
                           colors=[(0.5,0.5,0.5,1)]*2 + [color + (alpha,) for alpha in np.exp(-(np.logspace(0.20,2)-1.7)**2/(2.5**2*2.))], #smooth=3,
                           filled=True,
                           #linewidths=[1.0]*5,
                           layer='mask',
                           zorder=10, convention='calabretta')
            F.add_colorbar()
            F.colorbar.set_axis_label_text('T (K)')
            F.colorbar.set_axis_label_font(size=18)
            F.colorbar.set_label_properties(size=16)
            F.show_markers(sgrb2x, sgrb2y, color='k', facecolor='k', s=250,
                           edgecolor='k', alpha=0.9)
            F.save(os.path.join(figurepath, "big_maps", 'lores{0}{1}{2}_tmap_withmask.pdf'.format(smooth, outtype, vmax_str)))
            F.recenter(**big_recen)
            F.save(os.path.join(figurepath, "big_maps", 'big_lores{0}{1}{2}_tmap_withmask.pdf'.format(smooth, outtype, vmax_str)))
            log.info(os.path.join(figurepath, "big_maps", 'big_lores{0}{1}{2}_tmap_withmask.pdf'.format(smooth, outtype, vmax_str)))

            F.show_contour(dustcolumn,
                           levels=[5], colors=[(0,0,0,0.5)], zorder=15,
                           alpha=0.5,
                           linewidths=[0.5],
                           layer='dustcontour')
            F.recenter(**small_recen)
            F.save(os.path.join(figurepath, "big_maps", 'lores{0}{1}{2}_tmap_withcontours.pdf'.format(smooth, outtype, vmax_str)))
            F.recenter(**big_recen)
            F.save(os.path.join(figurepath, "big_maps", 'big_lores{0}{1}{2}_tmap_withcontours.pdf'.format(smooth, outtype, vmax_str)))
            log.info(os.path.join(figurepath, "big_maps", 'big_lores{0}{1}{2}_tmap_withcontours.pdf'.format(smooth, outtype, vmax_str)))

            F.hide_layer('mask')
            F.recenter(**small_recen)
            F.save(os.path.join(figurepath, "big_maps", 'lores{0}{1}{2}_tmap_nomask_withcontours.pdf'.format(smooth, outtype, vmax_str)))
            F.recenter(**big_recen)
            F.save(os.path.join(figurepath, "big_maps", 'big_lores{0}{1}{2}_tmap_nomask_withcontours.pdf'.format(smooth, outtype, vmax_str)))


            fig7 = pl.figure(7, figsize=figsize)
            fig7.clf()
            Fsn = aplpy.FITSFigure(peaksn, convention='calabretta', figure=fig7)
            Fsn.show_grayscale(vmin=0, vmax=10, stretch='linear', invert=True)
            Fsn.add_colorbar()
            Fsn.colorbar.set_axis_label_text('Peak S/N')
            Fsn.colorbar.set_axis_label_font(size=18)
            Fsn.colorbar.set_label_properties(size=16)
            Fsn.set_tick_labels_format('d.dd','d.dd')
            Fsn.recenter(**big_recen)
            Fsn.save(os.path.join(figurepath, "big_maps", 'big_lores{0}{1}{2}_peaksn.pdf'.format(smooth, outtype, vmax_str)))


            F.hide_layer('dustcontour')
            dusttemperature = '/Users/adam/work/gc/gcmosaic_temp_conv36.fits'
            F.show_contour(dusttemperature,
                           levels=[20,25],
                           colors=[(0,0,x,0.5) for x in [0.9,0.7,0.6,0.2]], zorder=20)
            F.recenter(**small_recen)
            F.save(os.path.join(figurepath, "big_maps",'lores{0}{1}{2}_tmap_withtdustcontours.pdf'.format(smooth, outtype, vmax_str)))
            F.recenter(**big_recen)
            F.save(os.path.join(figurepath, "big_maps",'big_lores{0}{1}{2}_tmap_withtdustcontours.pdf'.format(smooth, outtype, vmax_str)))
            log.info(os.path.join(figurepath, "big_maps",'big_lores{0}{1}{2}_tmap_withtdustcontours.pdf'.format(smooth, outtype, vmax_str)))


            im = fits.getdata(h2copath+ftemplate.format(smooth))
            data = im[np.isfinite(im)]

            fig9 = pl.figure(9)
            fig9.clf()
            ax9 = fig9.gca()
            h,l,p = ax9.hist(data, bins=np.linspace(0,300), alpha=0.5)
            shape, loc, scale = ss.lognorm.fit(data, floc=0)
            # from http://nbviewer.ipython.org/url/xweb.geos.ed.ac.uk/~jsteven5/blog/lognormal_distributions.ipynb
            mu = np.log(scale) # Mean of log(X) [but I want mean(x)]
            sigma = shape # Standard deviation of log(X)
            M = np.exp(mu) # Geometric mean == median
            s = np.exp(sigma) # Geometric standard deviation
            lnf = ss.lognorm(s=shape, loc=loc, scale=scale)
            pdf = lnf.pdf(np.arange(300))
            label1 = ("$\sigma_{{\mathrm{{ln}} x}} = {0:0.2f}$\n"
                      "$\mu_x = {1:0.2f}$\n"
                      "$\sigma_x = {2:0.2f}$".format(sigma, scale,s))
            pm = np.abs(ss.lognorm.interval(0.683, s=shape, loc=0, scale=scale) - scale)
            label2 = ("$x = {0:0.1f}^{{+{1:0.1f}}}_{{-{2:0.1f}}}$\n"
                      "$\sigma_{{\mathrm{{ln}} x}} = {3:0.1f}$\n"
                      .format(scale,
                              pm[1],
                              pm[0],
                              sigma,
                             ))
            ax9.plot(np.arange(300), pdf*h.max()/pdf.max(), linewidth=4, alpha=0.5,
                     label=label2)
            ax9.legend(loc='best')
            ax9.set_xlim(0,300)

            fig9.savefig(os.path.join(figurepath, "big_maps",
                                      'histogram_{0}{1}{2}_tmap.pdf'.format(smooth,
                                                                         outtype, vmax_str)),
                         bbox_inches='tight')

        #F.show_contour('h2co218222_all.fits', levels=[1,7,11,20,38], colors=['g']*5, smooth=1, zorder=5)
        #F.show_contour(datapath+'APEX_H2CO_merge_high_smooth_noise.fits', levels=[0.05,0.1], colors=['#0000FF']*2, zorder=3, convention='calabretta')
        #F.show_contour(datapath+'APEX_H2CO_merge_high_nhits.fits', levels=[9], colors=['#0000FF']*2, zorder=3, convention='calabretta',smooth=3)
        #F.show_regions('2014_expansion_targets_simpler.reg')
        #F.save('CMZ_H2CO_observed_planned.pdf')
        #F.show_rgb(background, wcs=wcs)
        #F.save('CMZ_H2CO_observed_planned_colorful.pdf')


fig = pl.figure(5, figsize=figsize)
fig.clf()
F2 = aplpy.FITSFigure(dusttemperature, convention='calabretta', figure=fig)
F2.show_colorscale(cmap=pl.cm.hot, vmin=10, vmax=40)
F2.add_colorbar()
F2.show_contour(h2copath+'H2CO_321220_to_303202_smooth_bl_integ_temperature.fits',
                convention='calabretta',
                levels=[30,75,100,150],
                cmap=pl.cm.BuGn)
F2.recenter(**small_recen)

F2.show_markers(sgrb2x, sgrb2y, color='k', facecolor='k', s=250,
               edgecolor='k', alpha=0.9)

F2.save(os.path.join(figurepath, "big_maps",'H2COtemperatureOnDust.pdf'))

F2.recenter(**big_recen)
F2.save(os.path.join(figurepath, "big_maps",'big_H2COtemperatureOnDust.pdf'))

fig = pl.figure(6, figsize=figsize)
fig.clf()
F = aplpy.FITSFigure('/Users/adam/work/gc/Tkin-GC.fits.gz',
                     convention='calabretta',
                     figure=fig)

for vmax in (100,200):
    cm = copy.copy(cmap)
    cm.set_bad((0.5,)*3)
    F.show_colorscale(cmap=cm,vmin=vmin,vmax=vmax)
    F.set_tick_labels_format('d.dd','d.dd')
    F.recenter(**small_recen)
    F.add_colorbar()
    F.colorbar.set_axis_label_text('T (K)')
    F.colorbar.set_axis_label_font(size=18)
    F.colorbar.set_label_properties(size=16)
    F.show_markers(sgrb2x, sgrb2y, color='k', facecolor='k', s=250,
                   edgecolor='k', alpha=0.9)

    F.save(os.path.join(figurepath, "big_maps", 'ott2014_nh3_tmap_15to{0}.pdf'.format(vmax)))
    F.show_colorscale(cmap=cm,vmin=vmin,vmax=80)
    F.save(os.path.join(figurepath, "big_maps", 'ott2014_nh3_tmap_15to80.pdf'))

    F.show_contour(dustcolumn,
                   levels=[5], colors=[(0,0,0,0.5)], zorder=15,
                   alpha=0.5,
                   linewidths=[0.5],
                   layer='dustcontour')
    F.save(os.path.join(figurepath, "big_maps", 'ott2014_nh3_tmap_15to80_withcontours.pdf'))
    F.show_colorscale(cmap=cm,vmin=vmin,vmax=vmax)
    F.save(os.path.join(figurepath, "big_maps", 'ott2014_nh3_tmap_15to{0}_withcontours.pdf'.format(vmax)))
