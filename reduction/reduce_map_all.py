"""

"""
from shfi_otf_pipeline import make_apex_cubes
from sdpy import plait
import FITS_tools
import os
from astropy.io import fits
from astropy import log
outdir = '/Users/adam/work/h2co/apex/april2014/'

def reduce_all_cubes_for_map(mapname, lowhigh='high', **kwargs):
    for dataset,maps in make_apex_cubes.datasets_2014.items():
        if mapname in maps:
            date = dataset[-10:]
            both_directions = True
            try:
                make_apex_cubes.build_cube_2014(mapname,
                                                lowhigh=lowhigh,
                                                posang=[50,70],
                                                datasets=[dataset],
                                                extra_suffix='_cal{0}_lscans'.format(date),
                                                **kwargs)
            except IndexError:
                both_directions = False

            try:
                make_apex_cubes.build_cube_2014(mapname,
                                                lowhigh=lowhigh,
                                                posang=[140,160],
                                                datasets=[dataset],
                                                extra_suffix='_cal{0}_bscans'.format(date),
                                                **kwargs)
            except IndexError:
                both_directions = False

            if both_directions:
                fileb = os.path.join(outdir, 'APEX_H2CO_2014_{1}_{2}_cal{0}_bscans.fits'.format(date, mapname, lowhigh))
                filel = os.path.join(outdir, 'APEX_H2CO_2014_{1}_{2}_cal{0}_lscans.fits'.format(date, mapname, lowhigh))

                cubeb = fits.getdata(fileb)
                cubel = fits.getdata(filel)

                if cubeb.shape != cubel.shape:
                    header = FITS_tools.fits_overlap(fileb, filel)
                    hdb = fits.getheader(fileb)
                    # Add back 3rd dimension... HACK
                    for key in hdb:
                        if key[0] == 'C' and key.strip()[-1] == '3':
                            header[key] = hdb[key]

                    FITS_tools.regrid_fits_cube(fileb, outheader=header, outfilename=fileb, clobber=True)
                    FITS_tools.regrid_fits_cube(filel, outheader=header, outfilename=filel, clobber=True)

                cubeb = fits.getdata(fileb)
                cubel = fits.getdata(filel)

                if cubeb.shape != cubel.shape:
                    log.fatal("Cube shapes don't match: {0}, {1}".format(cubeb.shape,cubel.shape))
                    raise ValueError

                cube_comb = plait.plait_cube([cubeb,cubel], angles=[0, 90], scale=5)
                cube_comb_naive = (cubeb+cubel)/2.

                header = fits.getheader(fileb)
                fits.PrimaryHDU(data=cube_comb,
                                header=header).writeto(os.path.join(outdir, '{1}_{2}_cal{0}_plait.fits'.format(date, mapname, lowhigh)),
                                                       clobber=True)
        #fits.PrimaryHDU(data=cube_comb_naive, header=header).writeto(os.path.join(outdir, 'MAP_001_high_cal{0}_naive.fits'), clobber=True)
        #fits.PrimaryHDU(data=cube_comb_naive-cube_comb, header=header).writeto(os.path.join(outdir, 'MAP_001_high_cal{0}_diff.fits'), clobber=True)
