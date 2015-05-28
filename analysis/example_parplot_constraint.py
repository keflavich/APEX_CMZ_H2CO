import pprint
import pylab as pl
import numpy as np
from constrain_parameters import paraH2COmodel
import paths

mf = paraH2COmodel()
mf.set_constraints(ratio321303=0.35, eratio321303=0.01, logh2column=23,
                   elogh2column=1, logabundance=np.log10(1.2e-9),
                   elogabundance=1, mindens=4, linewidth=10, taline303=0.571, etaline303=0.07,
                   taline321=0.2, etaline321=0.07)
pprint.pprint(mf.get_parconstraints(), width=1)


pl.figure(1)
mf.parplot1d('dens')
pl.figure(2)
mf.parplot1d('tem')
pl.figure(3)
mf.parplot1d('col')
pl.figure(4)
mf.parplot1d_all(levels=[0.68268949213708585])
pl.savefig(paths.fpath('param_fits/oned_fit_parameters_example.png'),
           bbox_inches='tight')
pl.figure(5)
mf.denstemplot()
pl.figure(7)
mf.denscolplot()
pl.figure(8)
mf.coltemplot()

pl.draw()
pl.show()
