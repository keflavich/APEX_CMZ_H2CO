"""
Determine what the real errors in the model are given errors on the parameters
by finding the best fit model for a series of measurements within the measured
error distribution
"""
import numpy as np
from constrain_parameters import paraH2COmodel
from astropy.utils.console import ProgressBar

model = paraH2COmodel()

ratio303321 = 0.3
eratio303321 = 0.01
logabundance = np.log10(1.2e-9)
elogabundance = 1.0
logh2column = 22.
elogh2column = 1.
mindens = 4.5
emindens = 0.2
linewidth=5.0

model.set_constraints(ratio303321=ratio303321,
                      eratio303321=eratio303321,
                      logabundance=logabundance,
                      elogabundance=elogabundance,
                      logh2column=logh2column,
                      elogh2column=elogh2column,
                      mindens=mindens,
                      emindens=emindens,
                      linewidth=linewidth)
default_constraints = model.get_parconstraints()

nsamples = 15000
random_samples = np.random.randn(nsamples)*eratio303321 + ratio303321

varnames = ('temperature','column','density')
results = {vn:[] for vn in varnames}
for ratio in ProgressBar(random_samples):
    model.chi2_r303321 = model.grid_getmatch_321to303(ratio,
                                                      eratio303321)

    model.compute_chi2_fromcomponents()
    constraints = model.get_parconstraints()

    for vn in varnames:
        results[vn].append(constraints['{0}_chi2'.format(vn)])

import pylab as pl
fig1 = pl.figure(1)
fig1.clf()
#ax1 = pl.subplot(311)
ax1 = fig1.gca()
p,l,h = ax1.hist(results['temperature'], histtype='stepfilled', color='r', bins=30)
ax1.vlines([default_constraints['tmin1sig_chi2'],
            default_constraints['tmax1sig_chi2']],
           0, np.max(p), linestyle='--', color='k')
#ax2 = pl.subplot(312)
#ax2.hist(results['column'], histtype='stepfilled')
#ax3 = pl.subplot(313)
#ax3.hist(results['density'], histtype='stepfilled')

fig2 = pl.figure(2)
fig2.clf()
model.parplot('dens','tem')

pl.draw()
pl.show()
