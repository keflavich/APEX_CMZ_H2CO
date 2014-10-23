"""
October 3, 2014:
Play with 'plait' method for cleaning up scan artifacts for a cube region.  Works great for MAP_003!

The 'scale' is actually a fourier-space scale, so smaller -> larger spatial scales

Next challenge: apply to 2013 and earlier data
"""
from shfi_otf_pipeline import make_apex_cubes
from sdpy import plait
import os
from astropy.io import fits
outdir = '/Users/adam/work/h2co/apex/april2014/'

make_apex_cubes.build_cube_2014('MAP_003',
                                lowhigh='high',
                                posang=[50,70],
                                datasets=[x
                                          for x,y in make_apex_cubes.datasets_2014.items()
                                          if 'MAP_003' in y],
                                extra_suffix='_lscans')

make_apex_cubes.build_cube_2014('MAP_003',
                                lowhigh='high',
                                posang=[140,160],
                                datasets=[x
                                          for x,y in make_apex_cubes.datasets_2014.items()
                                          if 'MAP_003' in y],
                                extra_suffix='_bscans')

fileb = os.path.join(outdir, 'APEX_H2CO_2014_MAP_003_high_bscans.fits')
filel = os.path.join(outdir, 'APEX_H2CO_2014_MAP_003_high_lscans.fits')

cubeb = fits.getdata(fileb)
cubel = fits.getdata(filel)
assert cubeb.shape == cubel.shape

cube_comb = plait.plait_cube([cubeb,cubel], angles=[0, 90], scale=5)
cube_comb_naive = (cubeb+cubel)/2.

header = fits.getheader(fileb)
fits.PrimaryHDU(data=cube_comb, header=header).writeto(os.path.join('MAP_003_high_plait.fits'), clobber=True)
fits.PrimaryHDU(data=cube_comb_naive, header=header).writeto(os.path.join(outdir, 'MAP_003_high_naive.fits'), clobber=True)
fits.PrimaryHDU(data=cube_comb_naive-cube_comb, header=header).writeto(os.path.join(outdir, 'MAP_003_high_diff.fits'), clobber=True)

import pylab as pl
pl.figure(1, figsize=(12,12))
pl.clf()
for ii,scale in enumerate((1,5,10)):
    cube_comb = plait.plait_cube([cubeb,cubel], angles=[0, 90], scale=scale)
    integ_plait = cube_comb[900:990,:,:].mean(axis=0)
    integ_naive = cube_comb_naive[900:990,:,:].mean(axis=0)
    integ_plait = cube_comb[2515:2605,:,:].mean(axis=0)
    integ_naive = cube_comb_naive[2515:2605,:,:].mean(axis=0)
    integ_diff = integ_plait-integ_naive
    pl.subplot(3,3,1+3*ii)
    pl.imshow(integ_plait)
    pl.title('Plait: Scale {0} pix'.format(scale))
    pl.subplot(3,3,2+3*ii)
    pl.imshow(integ_naive)
    pl.title('Naive')
    pl.subplot(3,3,3+3*ii)
    pl.imshow(integ_diff)
    pl.title('Diff: Plait-Naive')
    pl.savefig(os.path.join(outdir, 'PlaitVsNot_Comparison.pdf'), bbox_inches='tight')
