"""
Jan 15: This isn't tracked - why?  Is it really not needed?  I'm tracking it
now, but I think it is entirely redundant with make_ratio_cubes_ims and that
file should supercede this now.  BUT, make_ratio_cubes_ims doesn't include
ph2cogrid!  I guess BOTH are necessary.  make_ratio_cubes_ims is used for 321
cubes (which were a failed experiment from early December) and this is the
"real" one...


Make ratio maps and temperature maps of cubes and images (not cube sims)

try make_ratio_integ instead?
"""
from __future__ import print_function
import os
import numpy as np
from astropy.io import fits
from astropy import log

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
        mask = f[0].data.astype('bool')

        # mask now!
        f[0].data = ratio * mask

        f.writeto(h2copath+'H2CO_322221_to_303202_cube{0}.fits'.format(smooth),clobber=True)

        weight = top
        weight[weight < 0] = 0
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
        ratio_weighted = np.nansum(weight*ratio, axis=0) / np.nansum(weight, axis=0)
        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask_integ.fits'.format(smooth))
        f[0].data = ratio_weighted
        f.writeto(h2copath+'H2CO_321220_to_303202{0}_integ_weighted.fits'.format(smooth),clobber=True)
        ratio_masked_weighted = np.nansum(mask*weight*ratio, axis=0) / np.nansum(mask*weight, axis=0)
        f[0].data = ratio_masked_weighted
        f.writeto(h2copath+'H2CO_321220_to_303202{0}_integ_masked_weighted.fits'.format(smooth),clobber=True)

        log.info("Completed cube ratios using {0}".format(smooth))

def ph2cogrid(ntemp=50, trange=[10,200], abundances=(10**-8.5,10**-9),
              Nh2=(3e22,3e23), logdensities=(4,5)):
    import pyradex

    temperatures=np.linspace(trange[0],trange[1],ntemp)

    # initial density; will be modified later
    density = 1e4

    deltav = 5.0 # km/s

    R = pyradex.Radex(species='ph2co-h2',
                      abundance=abundances[0],
                      collider_densities={'H2':density},
                      deltav=deltav,
                      column=None,#Nh2[0]/abundances[0],
                      temperature=temperatures[0],
                      )

    Xarr = {}
    for abundance in abundances:
        Xarr[abundance] = {}
        for h2column in Nh2:

            densities = [10**x for x in logdensities]
            ratio1 = {d:[] for d in densities}
            ratio2 = {d:[] for d in densities}
            f1 = {d:[] for d in densities}
            f2 = {d:[] for d in densities}
            f3 = {d:[] for d in densities}

            for density in densities:
                R.density = {'H2': density}
                for temperature in temperatures:
                    R.abundance = abundance
                    R.temperature = temperature
                    log.info(str((R.run_radex(validate_colliders=False, reload_molfile=False),
                                  R.column, R.density['H2'],
                                  R.temperature, R.abundance)))

                    F1 = R.T_B[2]  # 218.222192 3_0_3
                    F2 = R.T_B[12] # 218.760066 3_2_1
                    F3 = R.T_B[9]  # 218.475632 3_2_2

                    ratio1[density].append(F2/F1)
                    ratio2[density].append(F3/F1)
                    f3[density].append(F3)
                    f2[density].append(F2)
                    f1[density].append(F1)
                print()

            f1 = {d:np.array([x.value for x in f1[d]]) for d in densities}
            f2 = {d:np.array([x.value for x in f2[d]]) for d in densities}
            f3 = {d:np.array([x.value for x in f3[d]]) for d in densities}
            ratio1 = {d:np.array(ratio1[d]) for d in densities}
            ratio2 = {d:np.array(ratio2[d]) for d in densities}

            Xarr[abundance][h2column] = {'flux1':f1,
                                         'flux2':f2,
                                         'flux3':f3,
                                         'ratio1':ratio1,
                                         'ratio2':ratio2}

    return Xarr

class TemperatureMapper(object):
    """
    For lazier evaluation of temperature mapping function
    """
    def __init__(self, trange=[10,300], ntemp=100, **kwargs):
        self.trange = trange
        self.ntemp = ntemp
        self.kwargs = kwargs

    def init(self):
        self.Xarr = ph2cogrid(trange=self.trange, ntemp=self.ntemp,
                              logdensities=(4,4.5,5), abundances=(1.2e-9,),
                              Nh2=(5e22,5e22), **self.kwargs)
        self.temperatures = np.linspace(self.trange[0], self.trange[1],
                                        self.ntemp)


    def get_mapper(self, lineid, tmin=np.nan, tmax=np.nan,
                   density=1e4,
                   column=5e22):
        if not hasattr(self,'temperatures'):
            self.init()

        rationame = {'321220': 'ratio1',
                     '322221': 'ratio2'}[lineid]

        # ugly hack because ph2co is indexed with floats
        # Use FIXED abundance, FIXED column, FIXED density
        ratios = self.Xarr[self.Xarr.keys()[0]][column][rationame][density]

        def ratio_to_tem(r):
            inds = np.argsort(ratios)
            return np.interp(r, ratios[inds], self.temperatures[inds], tmin,
                             tmax)

        return ratio_to_tem

    def __call__(self, x, lineid='321220', **kwargs):
        return self.get_mapper(lineid, **kwargs)(x)

if 'tm' not in locals():
    tm = TemperatureMapper(trange=[10,300],ntemp=100)

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
    temperaturemap(tm, path=h2copath, ratio=False, Nnsuffix='_dens1e5_col5e22',
                   density=1e5)
    temperaturemap(tm, path=h2copath, ratio=False, Nnsuffix='_dens3e4_col5e22',
                   density=10**4.5)
    temperaturemap(tm, path=h2copath, ratio=False, Nnsuffix='_dens1e4_col5e22',
                   density=1e4)
    # Higher column, higher density for Brick, 20/50 kms cloud, Sgr B2 regino
    temperaturemap(tm, path=h2copath, ratio=False, Nnsuffix='_dens1e5_col3e23',
                   density=1e5)

def temperaturemap(ratio_to_tem, path=h2copath, Nnsuffix="", ratio=True,
                   **kwargs):
    if ratio:
        doratio()

    import scipy.stats

    for suf_ in ('{0}','{0}_integ','{0}_integ_weighted'):#,'_cube{0}'):
        for smooth in ('','_smooth','_bl','_smooth_bl'):

            suf = suf_.format(smooth)

            for highline in ('321220','322221'):
                pfx = '{2}/H2CO_{0}_to_303202{1}'.format(highline,suf,path)
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
