import numpy as np
import os
import paths
import glob
import itertools
from astropy.io import fits
import FITS_tools
import scipy.stats
from sdpy import plait
from reduce_map_all import reduce_all_cubes_for_map
from astropy import log
from collections import defaultdict

for lowhigh in ('low','high'):
    reduce_all_cubes_for_map('MAP_006', lowhigh=lowhigh, calibration_factors=defaultdict(lambda: 1))

    files = glob.glob(os.path.join(paths.april2014path, '*MAP_006*scans.fits'))

    for fn1,fn2 in itertools.combinations(files, 2):
        header = FITS_tools.fits_overlap(fn1, fn2)

    hd = fits.getheader(fn1)
    # Add back 3rd dimension... HACK
    for key in hd:
        if key[0] == 'C' and key.strip()[-1] == '3':
            header[key] = hd[key]

    for fn in files:
        FITS_tools.regrid_fits_cube(fn, outheader=header, outfilename=fn, clobber=True)


    cube530b = fits.getdata(os.path.join(paths.april2014path, 'APEX_H2CO_2014_MAP_006_{0}_cal2014-05-30_bscans.fits'.format(lowhigh)))
    cube530l = fits.getdata(os.path.join(paths.april2014path, 'APEX_H2CO_2014_MAP_006_{0}_cal2014-05-30_lscans.fits'.format(lowhigh)))
    cube730b = fits.getdata(os.path.join(paths.april2014path, 'APEX_H2CO_2014_MAP_006_{0}_cal2014-07-31_bscans.fits'.format(lowhigh)))
    cube730l = fits.getdata(os.path.join(paths.april2014path, 'APEX_H2CO_2014_MAP_006_{0}_cal2014-07-30_lscans.fits'.format(lowhigh)))

    cube730comb = plait.plait_cube([cube730b,cube730l], angles=[0, 90], scale=5)
    cube530comb = plait.plait_cube([cube530b,cube530l], angles=[0, 90], scale=5)


    import pylab as pl
    fig2 = pl.figure(2)
    fig2.clf()
    ax2 = fig2.gca()
    ax2.grid()

    data = []
    for threshold in np.logspace(-0.6,1.1,10):
        mask = (cube530comb>threshold)&(cube730comb>threshold)
        if not np.any(mask):
            log.info("Skip threshold {0} for {1} because no data".format(threshold, lowhigh))
            continue
        ratio = (cube530comb/cube730comb)
        center,width = scipy.stats.norm.fit(ratio[mask])
        ax2.hist(ratio[mask],
                 bins=np.linspace(0,2,100), histtype='step', log=True, linewidth=2, alpha=0.7, label="{0:0.1f}".format(threshold))
        data.append([threshold,center,width])
    data = np.array(data)
    print(data)

    pl.figure(1)
    pl.clf()
    pl.grid()
    pl.errorbar(data[:,0], data[:,1], yerr=data[:,2], alpha=0.5, linewidth=2)
    pl.xlabel("Threshold")
    pl.ylabel("Ratio 05-30 / 07-30")
    pl.savefig(paths.fpath('calibration/MAP_006_0530vs0730_{0}_ratiovsthreshold.png'.format(lowhigh)), bbox_inches='tight')

    ax2.set_xlabel("Ratio 05-30 / 07-30")
    ax2.legend(loc='best')
    fig2.savefig(paths.fpath('calibration/MAP_006_0530vs0730_{0}_ratiohistograms.png'.format(lowhigh)), bbox_inches='tight')
