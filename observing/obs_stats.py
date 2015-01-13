from paths import rpath,mpath,opath,fpath
from astropy.io import fits
import pyregion
import numpy as np
from astropy import table
import paths
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

noise_img = fits.open(mpath('APEX_H2CO_merge_high_plait_all_noise.fits'))[0]
nhits_img = fits.open(mpath('APEX_H2CO_merge_high_plait_all_nhits.fits'))[0]

map_regions = pyregion.open(rpath('target_fields_8x8.reg'))

nhits = {}
noise = {}

for reg in map_regions:
    regs = pyregion.ShapeList([reg])
    mask = (regs.get_mask(hdu=noise_img) & np.isfinite(noise_img.data) &
            np.isfinite(nhits_img.data))

    name = reg.attr[1]['text'].split()[0]

    noise[name] = {'mean':noise_img.data[mask].mean(),
                   'median':np.median(noise_img.data[mask]),
                   'std':noise_img.data[mask].std()}

    nhits[name] = {'mean':nhits_img.data[mask].mean(),
                   'median':np.median(nhits_img.data[mask]),
                   'std':nhits_img.data[mask].std()}

keys = sorted(nhits.keys())

tbl = table.Table()
tbl.add_column(table.Column(data=keys, name='FieldID', dtype=str))
for dt,dtn in zip((noise,nhits), ("noise","nhits")):
    for ii,coltype in enumerate(('mean','median','std')):
        col = table.Column(data=[dt[name][coltype] for name in keys],
                           name=dtn+"_"+coltype, dtype='float')
        tbl.add_column(col)

tbl.write(opath('field_noise_stats.ipac'), format='ascii.ipac')
tbl.sort('noise_mean')
tbl[::-1].write(opath('field_noise_stats_sorted.ipac'), format='ascii.ipac')

import pylab as pl
pl.figure(1)
pl.clf()
ax = pl.gca()
ax.plot(tbl['nhits_mean'],tbl['noise_mean'], '.', zorder=5, markersize=15, alpha=0.75)
exptime = np.arange(*ax.get_xlim())
exptime = np.arange(60,220,dtype='float')
# Not clear if 0.75 is physical or a fit or what...
#pl.plot(exptime, 0.75/np.sqrt(np.arange(20,220,dtype='float')),
#        label='$1/\\sqrt{t}$', linewidth=2, alpha=0.5, color='k', zorder=-5)
# Apparently each "weight point" is 1/8s; this is empirical though
pl.plot(exptime, 2**0.5*155/(0.733e6 * exptime/8.)**0.5,
        label=r'$\sim1/\sqrt{t}$', linewidth=3, alpha=0.25, color='k',
        zorder=-5)
pl.axis((60,220,0.045,0.09))
pl.xlabel("Sum of 1/variance * gaussian weights")
pl.xlabel(r"Weighted Exposure Time ($\sim s$)", labelpad=20)
pl.ylabel("Noise per pixel in a 1 km s$^{-1}$ bin (K)")
pl.legend(loc='best')
pl.savefig(opath("observing_stats.png"), bbox_inches='tight')
pl.savefig(fpath("observing_stats.pdf"), bbox_inches='tight')
#pl.errorbar(tbl['nhits_mean'],tbl['noise_mean'],yerr=tbl['noise_std'],linestyle='none',marker='s')
#pl.plot(tbl['nhits_mean'],tbl['nhits_median'],'.')
