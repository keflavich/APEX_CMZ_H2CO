from astropy.io import fits
from matplotlib import cm
import aplpy
import paths

# create "coverage" mask map
covgf = fits.open(paths.mpath('APEX_H2CO_merge_high_nhits.fits'))
covgf[0].data = (covgf[0].data > 0.05).astype('int')
covgf.writeto(paths.mpath('coverage_bool.fits'),clobber=True)

F = aplpy.FITSFigure('/Users/adam/Dropbox/SMA_CMZ_FITS_files/gcmosaic_column_conv25.fits',convention='calabretta')
F.show_colorscale(cmap=cm.hot_r,vmax=100,vmin=0, vmid=-2.5, stretch='log')
F.recenter(0.6,-0.1,height=0.5,width=2.6)
F.show_contour(paths.mpath('coverage_bool.fits'),levels=[0.5], colors=['b'], convention='calabretta')
F.set_tick_labels_xformat('d.dd')
F.set_tick_labels_yformat('d.dd')
F.save(paths.fpath('coverage_on_herschel.png'),dpi=150)
F.save(paths.fpath('coverage_on_herschel.pdf'),dpi=150)
F.remove_layer('contour_set_1')
F.show_contour(paths.molpath('APEX_H2CO_303_202_smooth_bl_mask_integ.fits'),
               levels=[0.5,1,2], colors=['b'], convention='calabretta')
F.save(paths.fpath('H2CO_on_herschel.png'),dpi=150)
F.save(paths.fpath('H2CO_on_herschel.pdf'),dpi=150)
F.remove_layer('contour_set_2')
F.show_contour(paths.molpath('APEX_SiO_54_smooth_mask_integ.fits'),
               levels=[0.5,1,2], colors=['b'], convention='calabretta')
F.save(paths.fpath('SiO_on_herschel.png'),dpi=150)
F.save(paths.fpath('SiO_on_herschel.pdf'),dpi=150)
