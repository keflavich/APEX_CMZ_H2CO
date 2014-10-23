from paths import rpath,mpath,opath
from astropy.io import fits
import pyregion
import numpy as np
from astropy import table

noise_img = fits.open(mpath('APEX_H2CO_merge_high_plait_all_noise.fits'))[0]
nhits_img = fits.open(mpath('APEX_H2CO_merge_high_nhits.fits'))[0]

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
pl.plot(tbl['nhits_mean'],tbl['noise_mean'],'.', zorder=5)
# Not clear if 0.75 is physical or a fit or what...
pl.plot(np.arange(20,200), 0.75/np.sqrt(np.arange(20,200,dtype='float')),
        label='$1/\\sqrt{t}$', linewidth=2, alpha=0.5, color='k', zorder=-5)
pl.xlabel("Average number of 0.25s integrations per pixel")
pl.ylabel("Noise per pixel in a 1 km/s bin (K)")
pl.legend(loc='best')
pl.savefig(opath("observing_stats.png"))
#pl.errorbar(tbl['nhits_mean'],tbl['noise_mean'],yerr=tbl['noise_std'],linestyle='none',marker='s')
#pl.plot(tbl['nhits_mean'],tbl['nhits_median'],'.')
