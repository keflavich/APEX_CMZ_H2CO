import numpy as np
import os
from astropy.io import fits
from sdpy import plait
mergepath='/Volumes/passport/apex/merged_datasets'
mergefile2='APEX_H2CO_merge_high'

all_targets = ("_2014_bscans", "_2014_lscans", "_2013","_ao")
nplaittargets = 2
plait_targets = all_targets[:nplaittargets]
headers = [fits.getheader(os.path.join(mergepath, mergefile2+suff+".fits"))
           for suff in all_targets]
header = headers[0]
for h in headers:
    header.update(h)

cubes = [fits.getdata(os.path.join(mergepath, mergefile2+suff+".fits"))
         for suff in all_targets]
angles = [0, 90, 90-58.6, 90-58.6]
slices = [c[261:271,:,:].mean(axis=0) for c in cubes]
for s in slices:
    s[s==0] = np.nan

weights = [fits.getdata(os.path.join(mergepath, mergefile2+suff+"_nhits.fits"))
           for suff in all_targets]
sweights = [np.median(w[w>0]) for w in weights]

naive_wt = np.nansum([s*w for s,w in zip(slices,weights)[:]],
                      axis=0) / np.sum(weights[:],0)
naive = np.nanmean(slices, 0)

#cube_comb = plait.plait_cube(cubes, angles=angles, scale=3)
scale = 3
slice_comb = plait.plait_plane(slices[:nplaittargets], angles=angles[:nplaittargets],
                               scale=scale, weights=None)
slice_comb_wt = plait.plait_plane(slices[:nplaittargets],
                                  angles=angles[:nplaittargets], scale=scale,
                                  weights=sweights[:nplaittargets])
slice_comb_all = np.nansum([slice_comb*(weights[0]+weights[1]),
                            slices[2]*weights[2],
                            slices[3]*weights[3]],axis=0) / np.sum(weights,axis=0)
#naive[np.all(np.isnan(slices, 0), 0)] = np.nan
diff = naive_wt-slice_comb_all

#hdu = fits.PrimaryHDU(data=cube_comb, header=header)
#hdu.writeto(os.path.join(mergepath, mergefile2+"_plait.fits"), clobber=True)

import pylab as pl
pl.figure(1, figsize=(18,18))
pl.clf()
pl.subplot(5,1,1)
pl.imshow(slice_comb, vmin=-1, vmax=1)
pl.subplot(5,1,2)
pl.imshow(slice_comb_all, vmin=-1, vmax=1)
pl.subplot(5,1,3)
pl.imshow(diff, vmin=-1, vmax=1)
pl.subplot(5,1,4)
pl.imshow(naive, vmin=-1, vmax=1)
pl.subplot(5,1,5)
pl.imshow(naive_wt, vmin=-1, vmax=1)

pl.draw()
pl.show()
