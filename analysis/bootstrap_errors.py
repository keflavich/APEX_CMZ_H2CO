"""
Determine what the real errors in the model are given errors on the parameters
by finding the best fit model for a series of measurements within the measured
error distribution

EDIT 5/3/2015: After some further thinking, this approach doesn't really make sense.
We already have the likelihood of the model given the data.  This computes the
distribution of the maximum-likelihood model, which simply isn't the same
thing.

However, I think we *can* take the modeled parconstraints and sample from
within those to see what ratios we get out and whether that distribution matches
the input distribution.  
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
# assume the asymmetric errors correspond to symmetric gaussians in all cases
# (this is not a great assumption but is the simplest in context; really we're
# discussing some kind of ill-defined skewed distribution.  It's not a
# lognormal.)
rand_t = np.random.randn(nsamples)
err_t = ((default_constraints['temperature_chi2'] - default_constraints['tmin1sig_chi2']) +
         (default_constraints['tmax1sig_chi2'] - default_constraints['temperature_chi2']))/2.
temperatures = (default_constraints['temperature_chi2'] + (rand_t * err_t))

rand_c = np.random.randn(nsamples)
err_c = ((default_constraints['column_chi2'] - default_constraints['cmin1sig_chi2']) +
         (default_constraints['cmax1sig_chi2'] - default_constraints['column_chi2']))/2.
columns = (default_constraints['column_chi2'] + (rand_c * err_c))

rand_d = np.random.randn(nsamples)
err_d = ((default_constraints['density_chi2'] - default_constraints['dmin1sig_chi2']) +
         (default_constraints['dmax1sig_chi2'] - default_constraints['density_chi2']))/2.
densitys = (default_constraints['density_chi2'] + (rand_d * err_d))

tem_xarr = model.temparr[:,0,0]
den_xarr = model.densityarr[0,:,0]
col_xarr = model.columnarr[0,0,:]

ratios = []

for col, tem, den in zip(columns, temperatures, densitys):
    bestfit_pos = (np.argmin(abs(tem_xarr-tem)),
                   np.argmin(abs(den_xarr-den)),
                   np.argmin(abs(col_xarr-col)))
    ratio = model.modelratio1[bestfit_pos]
    ratios.append(ratio)


import pylab as pl
fig1 = pl.figure(1)
fig1.clf()
#ax1 = pl.subplot(311)
ax1 = fig1.gca()
p,l,h = ax1.hist(ratios, histtype='stepfilled', color='r', bins=60)
ax1.vlines([ratio303321-eratio303321,
            ratio303321,
            ratio303321+eratio303321],
           0, np.max(p), linestyle='--', color='k')
ax1.vlines([np.percentile(ratios, 15.8),
            np.percentile(ratios,100-15.8)],
           0, np.max(p), linestyle=':', color='k')

pl.draw()
pl.show()


# A different experiment: try sampling from the full likelihood array
# We need to convert the likelihoods to a CDF first...
from scipy import stats
# 3 "degrees of freedom"
cdf = stats.chi2.cdf(model.chi2 - model.chi2.min(), 3)
cdf1 = stats.chi2.cdf(model.chi2 - model.chi2.min(), 1)
#sortinds = np.argsort(cdf.ravel())
#sorted_cdf = cdf.flat[sortinds]

ratios = []

for rand in np.random.rand(nsamples):
    bestfit_pos = np.unravel_index(np.argmin(abs(cdf-rand)),
                                   model.chi2.shape)
    ratio = model.modelratio1[bestfit_pos]
    ratios.append(ratio)

fig2 = pl.figure(2)
fig2.clf()
#ax1 = pl.subplot(311)
ax2 = fig2.gca()
p,l,h = ax2.hist(ratios, histtype='stepfilled', color='r', bins=60)
ax2.vlines([ratio303321-eratio303321,
            ratio303321,
            ratio303321+eratio303321],
           0, np.max(p), linestyle='--', color='k')
ax2.vlines([np.percentile(ratios, 15.8),
            np.percentile(ratios, 50),
            np.percentile(ratios,100-15.8)],
           0, np.max(p), linestyle=':', color='k')

pl.draw()
pl.show()


"""
# This stuff is not entirely sensible: it's attempting to answer an ill-defined
# question.
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
"""
