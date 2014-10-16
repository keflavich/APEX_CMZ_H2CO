"""
Code to reduce cross-scans of a single field.  Written for testing, but can be
generally useful.
"""
from shfi_otf_pipeline import make_apex_cubes
from sdpy import plait
import os
from astropy.io import fits

outdir = make_apex_cubes.april2014path

def reduce_crosscans(map_name, lowhigh='high'):
    make_apex_cubes.build_cube_2014(map_name,
                                    lowhigh=lowhigh,
                                    posang=[50,70],
                                    datasets=[x
                                              for x,y in make_apex_cubes.datasets_2014.items()
                                              if map_name in y],
                                    extra_suffix='_lscans')

    make_apex_cubes.build_cube_2014(map_name,
                                    lowhigh=lowhigh,
                                    posang=[140,160],
                                    datasets=[x
                                              for x,y in make_apex_cubes.datasets_2014.items()
                                              if map_name in y],
                                    extra_suffix='_bscans')

    fileb = os.path.join(outdir,
                         'APEX_H2CO_2014_{0}_{1}_bscans.fits'.format(map_name,
                                                                     lowhigh))
    filel = os.path.join(outdir,
                         'APEX_H2CO_2014_{0}_{1}_lscans.fits'.format(map_name,
                                                                     lowhigh))

    cubeb = fits.getdata(fileb)
    cubel = fits.getdata(filel)
    assert cubeb.shape == cubel.shape

    cube_comb = plait.plait_cube([cubeb,cubel], angles=[0, 90], scale=3)

    header = fits.getheader(fileb)
    hdu = fits.PrimaryHDU(data=cube_comb, header=header)
    hdu.writeto(os.path.join(outdir, '{0}_{1}_plait.fits'.format(map_name, lowhigh)),
                clobber=True)
