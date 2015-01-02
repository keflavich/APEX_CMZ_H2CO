import numpy as np
import matplotlib
from astropy.io import fits
import aplpy
from pyspeckit_fit_region import (do_1comp_region, do_g047_box,
                                  do_pyspeck_fits_1comp, do_pyspeck_fits_2comp,
                                  do_the_brick, get_subregion_pcube,
                                  remove_bad_pars)
import piecewise_rtotem
from paths import (h2copath, mergepath, figurepath, regpath, analysispath,
                   mpath, dpath)
import os

def fpath(x, figurepath=os.path.join(figurepath, 'pyspeckitmaps')):
    return os.path.join(figurepath, x)
def hpath(x, figurepath=os.path.join(h2copath, 'pyspeckitmaps')):
    return os.path.join(figurepath, x)

if not os.path.exists(fpath("")):
    os.mkdir(fpath(""))
if not os.path.exists(hpath("")):
    os.mkdir(hpath(""))

cm = matplotlib.cm.RdYlBu_r
cm.set_bad('#888888')

regions = {'G1.12-0.10': {'vrange':[-40,0], 'startpoint': (5,5)},
           'G1.00-0.02box': dict(vrange=[60,100], startpoint=(5,5)),
           'G1.00-0.11': dict(vrange=[100,140], startpoint=(5,4)),
           'G0.67-0.10': dict(vrange=[00,50], startpoint=(5,5)),
           'G0.47-0.07box': {'vrange':[50,125], 'startpoint':(22,6)},
           'G1.32-0.02box': {'vrange':[75,125], 'startpoint':(10,6)},
           'G1.59+0.01box': {'vrange':[130,200], 'startpoint':(13,28),
                             'minpeak':0.15, 'signal_cut': 2},
           'G1.21+0.00box': {'vrange':[60,100], 'startpoint':(5,5)},
           'G1.16-0.14': {'vrange':[120,160], 'startpoint':(5,5)}, # PROBABLY ALL BAD
          }

def fit_and_save(regname):
    pc = do_1comp_region(regname, **regions[regname])
    pc.mapplot(estimator=3, cmap=cm)
    
    plane = pc.mapplot.plane
    plane[plane == 0.0] = np.nan
    temmap = piecewise_rtotem.pwtem(plane.flat).reshape(plane.shape)
    hdu = fits.PrimaryHDU(data=temmap,
                          header=pc.mapplot.header)
    hdu.writeto(hpath("{0}_pyspeckit_tmap.fits".format(regname)),
                clobber=True)

    pc.mapplot.figure.clf()
    FF = aplpy.FITSFigure(hdu, figure=pc.mapplot.figure,
                          convention='calabretta')
    FF.set_tick_labels_format('d.dd','d.dd')
    FF.show_colorscale(cmap=cm, vmin=15, vmax=200, interpolation='nearest')
    FF.add_colorbar()
    FF.save(fpath("{0}_pyspeckit_tmap.pdf".format(regname)), dpi=150)
    #pl.savefig(fpath("{0}_pyspeckit_tmap.pdf".format(regname)))

    return pc

def fit_all(regions):
    for regname in regions:
        fit_and_save(regname)

