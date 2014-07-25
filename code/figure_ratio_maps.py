from paths import h2copath,figurepath
import os
import aplpy
import pylab as pl
import matplotlib

cm = matplotlib.cm.RdYlBu_r
cm.set_bad('#888888')

for bl in ("","_bl"):
    for smooth in ("","_smooth","_vsmooth"):
        ratio1 = 'H2CO_321220_to_303202{0}{1}_integ.fits'.format(smooth,bl)
        ratio2 = 'H2CO_322221_to_303202{0}{1}_integ.fits'.format(smooth,bl)


        for ii,ratio in enumerate((ratio1, ratio2)):
            fig = pl.figure(ii+1)
            fig.clf()
            F = aplpy.FITSFigure(os.path.join(h2copath, ratio),
                                 convention='calabretta', figure=fig)
            F.show_colorscale(cmap=cm)
            F.add_colorbar()
            F.tick_labels.set_xformat('d.dd')
            F.tick_labels.set_yformat('d.dd')
            F.recenter(0.31, -0.05, width=1.1, height=0.3)
            F.save(os.path.join(figurepath, ratio.replace(".fits",".pdf")))
