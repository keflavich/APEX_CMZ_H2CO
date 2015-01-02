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
           'G1.00-0.02': dict(vrange=[20,130], startpoint=(5,5)),
           'G0.67-0.10': dict(vrange=[00,50], startpoint=(5,5)),
           'G0.47-0.07box': {'vrange':[50,125], 'startpoint':(22,6)},
          }

for regname in regions:
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
    FF.show_colorscale(cmap=cm, vmin=15, vmax=200, interpolation='nearest')
    FF.add_colorbar()
    FF.save(fpath("{0}_pyspeckit_tmap.pdf".format(regname)), dpi=150)
    #pl.savefig(fpath("{0}_pyspeckit_tmap.pdf".format(regname)))
