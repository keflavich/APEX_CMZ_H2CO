import numpy as np
import pyspeckit
import paths
import make_apex_cubes
import pylab
pylab.ion()

name = 'MAP_001'
window = 'AP-H201-X202'
date = '2014-04-02'
sp = pyspeckit.Spectrum(paths.dpath('field_spectra/{name}_{window}_{date}.fits'.format(name=name,
                                                                                  window=window,
                                                                                  date=date)))
sp2 = pyspeckit.Spectrum(paths.dpath('field_spectra/{name}_{window}_{date}_fflagged.fits'.format(name=name,
                                                                                            window=window,
                                                                                            date=date)))

diff = (sp2-sp)
diff.specname=''
diff.xarr.convert_to_unit('GHz')
diff.plotter(ymax=-0.17, ymin=-0.3, alpha=1, xmin=217.500, xmax=218.500, linestyle='none', marker='.', markersize=2)
diff.plotter.savefig(paths.fpath('worst_baselines_map001.pdf'))

sp.xarr.refX = make_apex_cubes.bright_lines['H2CO_303_202']
sp.xarr.refX_unit = 'GHz'
sp.xarr.convert_to_unit('km/s')
sp.data = np.exp(-(sp.xarr-0)**2/(2*5.**2)) * 0.07
sp.xarr.convert_to_unit('GHz')
sp.specname=''
sp.plotter(axis=diff.plotter.axis, clear=False, color='b', offset=-0.25,
           ymax=-0.17, ymin=-0.3, xmin=217.500, xmax=218.500, 
           alpha=0.4, linewidth=2, zorder=-5
          )
diff.plotter.axis.set_ylim(-0.3, -0.17)
pylab.draw()
diff.plotter.savefig(paths.fpath('worst_baselines_map001_withsynth.pdf'))
