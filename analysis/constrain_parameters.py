"""
Functions for fitting temperature (and density and column) from the line ratio
plus whatever other constraints are available
"""
import numpy as np
from astropy import units as u
from h2co_modeling import grid_fitter

class paraH2COmodel(object):

    def __init__(self, tbackground=2.73, gridsize=250):
        from pyspeckit_fitting import (texgrid303, taugrid303, texgrid321, taugrid321,
                                       texgrid322, taugrid322, hdr)

        self.texgrid303 = texgrid303
        self.taugrid303 = taugrid303
        self.texgrid321 = texgrid321
        self.taugrid321 = taugrid321
        self.texgrid322 = texgrid322
        self.taugrid322 = taugrid322
        self.hdr = hdr

        self.Tbackground = tbackground
        self.tline303a = ((1.0-np.exp(-np.array(self.taugrid303))) *
                          (self.texgrid303-self.Tbackground))
        self.tline321a = ((1.0-np.exp(-np.array(self.taugrid321))) *
                          (self.texgrid321-self.Tbackground))
        self.tline322a = ((1.0-np.exp(-np.array(self.taugrid322))) *
                          (self.texgrid322-self.Tbackground))

        zinds,yinds,xinds = np.indices(self.tline303a.shape)
        upsample_factor = np.array([250./self.tline303a.shape[0],
                                    250./self.tline303a.shape[1],
                                    250./self.tline303a.shape[2]], dtype='float')
        uzinds,uyinds,uxinds = upsinds = np.indices([x*us
                                                     for x,us in zip(self.tline303a.shape,
                                                                     upsample_factor)],
                                                   dtype='float')
        self.tline303 = map_coordinates(self.tline303a,
                                   upsinds/upsample_factor[:,None,None,None],
                                   mode='nearest')
        self.tline321 = map_coordinates(self.tline321a,
                                   upsinds/upsample_factor[:,None,None,None],
                                   mode='nearest')
        self.tline322 = map_coordinates(self.tline322a,
                                   upsinds/upsample_factor[:,None,None,None],
                                   mode='nearest')
    
        self.tline = {303: self.tline303,
                      321: self.tline321,
                      322: self.tline322}

        self.densityarr = ((uxinds + self.hdr['CRPIX1']-1)*self.hdr['CDELT1'] /
                      float(upsample_factor[2])+self.hdr['CRVAL1']) # log density
        self.columnarr  = ((uyinds + self.hdr['CRPIX2']-1)*self.hdr['CDELT2'] /
                      float(upsample_factor[1])+self.hdr['CRVAL2']) # log column
        self.temparr    = ((uzinds + self.hdr['CRPIX3']-1)*self.hdr['CDELT3'] /
                      float(upsample_factor[0])+self.hdr['CRVAL3']) # lin temperature
        self.drange = [self.densityarr.min(), self.densityarr.max()]
        self.crange = [self.columnarr.min(),  self.columnarr.max()]
        self.trange = [self.temparr.min(),    self.temparr.max()]
        self.darr = self.densityarr[0,0,:]
        self.carr = self.columnarr[0,:,0]
        self.tarr = self.temparr[:,0,0]

        # While the individual lines are subject to filling factor uncertainties, the
        # ratio is not.
        self.modelratio1 = self.tline321/self.tline303
        self.modelratio2 = self.tline322/self.tline321

        self.model_logabundance = np.log10(10**self.columnarr / u.pc.to(u.cm) /
                                           10**self.densityarr)

    def grid_getmatch_321to303(self, ratio, eratio):
            match,indbest,chi2r = grid_fitter.grid_getmatch(ratio, eratio,
                                                            self.modelratio1)
            return match, indbest, chi2r

    def grid_getmatch_322to321(self, ratio, eratio):
            match,indbest,chi2r = grid_fitter.grid_getmatch(ratio, eratio,
                                                            self.modelratio2)
            return match, indbest, chi2r

    def chi2_fillingfactor(self, tline, etline, lineid):
        """
        Return a chi^2 value for each model parameter treating the specified
        line brightness as a lower limit

        Parameters
        ----------
        tline : float
            The line brightness temperature
        lineid : int
            The line id, one of 303,321,322
        """
        chi2 = ((self.tline[lineid] - tline)/etline)**2 * (self.tline[lineid] < tline)
        return chi2

    def chi2_column(self, logh2column, elogh2column, h2coabundance, linewidth):

        h2fromh2co = np.log10(10**self.columnarr * (np.sqrt(np.pi) * linewidth)
                              / 10**h2coabundance)
        chi2_h2 = ((h2fromh2co-logh2column)/elogh2column)**2

        return chi2_h2

    def chi2_abundance(self, logabundance, elogabundance):
        model_logabundance = np.log10(10**self.columnarr / u.pc.to(u.cm) /
                                      10**self.densityarr)
        chi2X = ((model_logabundance-logabundance)/elogabundance)**2
        return chi2X
