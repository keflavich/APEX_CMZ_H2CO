"""
Make ratio maps and temperature maps of cubes and images (not cube sims)
"""
from __future__ import print_function
import os
import numpy as np
from astropy.io import fits
from astropy import log

from temperature_mapper import ph2cogrid, TemperatureMapper, tm
from paths import h2copath

def doratio(h2copath=h2copath, maxratio=1):
    """
    I swapped top and bottom because that's what the models were set up for...
    this means 303 is on bottom, which it should be because it's the stronger line

    maxratio : float
        0.64 is the highest physical value... not really though.  If they're
        optically thick, they will of course saturate at 1.0

    """
    ##### integrated #####

    for smooth in ('','_smooth','_bl','_smooth_bl'):
        top = fits.getdata(h2copath+'APEX_H2CO_303_202{0}_mask_integ.fits'.format(smooth))
        bottom = fits.getdata(h2copath+'APEX_H2CO_322_221{0}_CH3OHchomped_mask_integ.fits'.format(smooth))

        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>maxratio] = np.nan

        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask_integ.fits'.format(smooth))

        f[0].data = ratio

        f.writeto(h2copath+'H2CO_322221_to_303202{0}_integ.fits'.format(smooth),clobber=True)

        bottom = fits.getdata(h2copath+'APEX_H2CO_321_220{0}_mask_integ.fits'.format(smooth))
        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>maxratio] = np.nan

        f[0].data = ratio

        f.writeto(h2copath+'H2CO_321220_to_303202{0}_integ.fits'.format(smooth),clobber=True)
        log.info("Completed integrated ratio maps using {0}".format(smooth))



    ##### cube #####
    for smooth in ('','_smooth','_bl','_smooth_bl'):
        top = fits.getdata(h2copath+'APEX_H2CO_303_202{0}.fits'.format(smooth))
        bottom = fits.getdata(h2copath+'APEX_H2CO_322_221{0}_CH3OHchomped_masked.fits'.format(smooth))

        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>maxratio] = np.nan

        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask.fits'.format(smooth))
        cubeheader = f[0].header
        mask = f[0].data.astype('bool')

        # mask now!
        f[0].data = ratio * mask

        f.writeto(h2copath+'H2CO_322221_to_303202_cube{0}.fits'.format(smooth),clobber=True)

        weight = top
        weight[weight < 0] = 0
        weight[~mask] = 0
        ratio_weighted = np.nansum(weight*ratio, axis=0) / np.nansum(weight, axis=0)
        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask_integ.fits'.format(smooth))
        f[0].data = ratio_weighted
        f.writeto(h2copath+'H2CO_322221_to_303202{0}_integ_weighted.fits'.format(smooth),clobber=True)
        ratio_masked_weighted = np.nansum(mask*weight*ratio, axis=0) / np.nansum(mask*weight, axis=0)
        f[0].data = ratio_masked_weighted
        f.writeto(h2copath+'H2CO_322221_to_303202{0}_integ_masked_weighted.fits'.format(smooth),clobber=True)


        bottom = fits.getdata(h2copath+'APEX_H2CO_321_220{0}.fits'.format(smooth))

        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>maxratio] = np.nan

        #f = fits.open(h2copath+'APEX_H2CO_303_202{0}.fits'.format(smooth))
        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask.fits'.format(smooth))

        # mask out...
        f[0].data = ratio * mask

        f.writeto(h2copath+'H2CO_321220_to_303202_cube{}.fits'.format(smooth),clobber=True)

        weight = top
        weight[weight < 0] = 0
        weight[~mask] = 0
        ratio_weighted = np.nansum(weight*ratio, axis=0) / np.nansum(weight, axis=0)
        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask_integ.fits'.format(smooth))
        f[0].data = ratio_weighted
        f.writeto(h2copath+'H2CO_321220_to_303202{0}_integ_weighted.fits'.format(smooth),clobber=True)
        ratio_masked_weighted = np.nansum(mask*weight*ratio, axis=0) / np.nansum(mask*weight, axis=0)
        f[0].data = ratio_masked_weighted
        f.writeto(h2copath+'H2CO_321220_to_303202{0}_integ_masked_weighted.fits'.format(smooth),clobber=True)

        # PV diagrams
        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask_integ.fits'.format(smooth))
        for kw in ['CRVAL','CDELT','CUNIT','CTYPE','CRPIX']:
            f[0].header[kw+"2"] = cubeheader[kw+"3"]
        f[0].header['CTYPE1'] = ('OFFSET','Galactic Longitude') # intentionally NOT a valid ctype
        del f[0].header['LONPOLE']
        del f[0].header['LATPOLE']
         
        weight[weight < 0] = 0
        weight[~mask] = 0
        ratio_weighted = np.nansum(weight*ratio, axis=1) / np.nansum(weight, axis=1)
        f[0].data = ratio_weighted
        f.writeto(h2copath+'pv_H2CO_321220_to_303202{0}_integ_weighted.fits'.format(smooth),clobber=True)

        ratio_masked_weighted = np.nansum(mask*weight*ratio, axis=1) / np.nansum(mask*weight, axis=1)
        f[0].data = ratio_masked_weighted
        f.writeto(h2copath+'pv_H2CO_321220_to_303202{0}_integ_masked_weighted.fits'.format(smooth),clobber=True)

        log.info("Completed cube ratios using {0}".format(smooth))


def do_temperature(ratio=True, h2copath=h2copath):
    #temperaturemap(tm, path=h2copath, ratio=ratio) # Defaults - unlabeled, useless
    # Why 5e22?
    # import masked_cubes
    # from higal_gridded import column_regridded
    # flatmask = masked_cubes.bmask._mask.max(axis=0)
    # np.median(column_regridded.data[flatmask])
    #
    # Why 3 densities?  There is probably a range from 1e4-1e5
    #
    temperaturemap(tm, path=h2copath, ratio=False, Nnsuffix='_dens1e5',
                   density=1e5)
    temperaturemap(tm, path=h2copath, ratio=False, Nnsuffix='_dens3e4',
                   density=10**4.5)
    temperaturemap(tm, path=h2copath, ratio=False, Nnsuffix='_dens1e4',
                   density=1e4)

def temperaturemap(ratio_to_tem, path=h2copath, Nnsuffix="", ratio=True,
                   **kwargs):
    if ratio:
        doratio()

    import scipy.stats

    for pre_ in ("pv_",""):
        for suf_ in ('{0}','{0}_integ','{0}_integ_weighted',"{0}_integ_masked_weighted"):#,'_cube{0}'):
            for smooth in ('','_smooth','_bl','_smooth_bl'):

                suf = suf_.format(smooth)

                for highline in ('321220','322221'):
                    pfx = '{2}/{3}H2CO_{0}_to_303202{1}'.format(highline,suf,path,pre_)
                    if os.path.exists(pfx+'.fits'):
                        log.info("Temperature mapping {0}".format(pfx))
                    else:
                        log.info("Skipping {0}".format(pfx))
                        continue
                    rmap = fits.getdata(pfx+'.fits')
                    tmap = ratio_to_tem(rmap, highline, **kwargs)#, tmin=10, tmax=300)
                    rf = fits.open(pfx+'.fits')
                    rf[0].header['BUNIT'] = 'K'
                    rf[0].header['BTYPE'] = 'TKIN'
                    del rf[0].header['LONPOLE']
                    del rf[0].header['LATPOLE']
                    rf[0].data = tmap
                    rf.writeto(pfx+'_temperature{suffix}.fits'.format(suffix=Nnsuffix),clobber=True)
                    if pre_ == "":
                        if 'cube' not in suf:
                            #mask = fits.getdata(path+'APEX_H2CO_322_221{0}_mask_integ.fits'.format(smooth)) > 0.0025
                            mask = fits.getdata(path+'APEX_H2CO_303_202{0}_mask_integ.fits'.format(smooth)) > 1.25
                        else:
                            mask = fits.getdata(path+'APEX_H2CO_303_202{0}_mask.fits'.format(smooth.rstrip('_bl'))).astype('bool')
                        tmap[~mask] = np.nan
                        rf[0].data = tmap
                        rf.writeto(pfx+'_temperature{suffix}_masked.fits'.format(suffix=Nnsuffix),
                                   clobber=True)

                    if 'cube' in suf:
                        rf[0].data = np.nanmax(tmap, axis=0)
                        rf.writeto(pfx+'_peaktemperature{suffix}.fits'.format(suffix=Nnsuffix), clobber=True)
                        rf[0].data = scipy.stats.nanmedian(tmap, axis=0)
                        rf.writeto(pfx+'_midtemperature{suffix}.fits'.format(suffix=Nnsuffix), clobber=True)

                        wt = fits.getdata(path+'APEX_H2CO_303_202{0}.fits'.format(smooth)).astype('bool')
                        wt[(True-mask) | (tmap==0) | (True-np.isfinite(tmap))] = 0
                        rf[0].data = np.nansum(tmap*wt, axis=0)/np.nansum(wt, axis=0)
                        rf.writeto(pfx+'_wtdmeantemperature{suffix}.fits'.format(suffix=Nnsuffix), clobber=True)

if __name__ == '__main__':
    doratio()
    do_temperature(ratio=False)
