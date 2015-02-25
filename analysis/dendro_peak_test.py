"""
Example showing that the dendrogram-extracted line width (or size scale) is not
representative of the "true" width
"""
import numpy as np
import pylab as pl
import paths

width = 1.3
center = 1.5
x = np.linspace(-10,10,1000)
y = np.exp(-(x-center)**2/(2.*width**2))

thresholds = np.logspace(-2,-0.05,100)
variances = []
for threshold in thresholds:
    mask = y > (y.max()*threshold)
    m0 = y[mask].sum()
    m1 = (x*y)[mask].sum()/m0
    m2 = (y[mask]*(x[mask]-m1)**2).sum()/m0
    variances.append(m2)

pl.figure(1)
pl.clf()
pl.plot(thresholds, np.sqrt(variances), alpha=0.5, linewidth=2)
pl.axhline(width, linestyle='--', color='k', alpha=0.5)
pl.ylabel("Variance measured")
pl.xlabel("Threshold level (fraction of peak)")
pl.draw()
pl.show()
pl.savefig(paths.fpath("dendro/dendro_leaf_variance_bias.pdf"), bbox_inches='tight')
