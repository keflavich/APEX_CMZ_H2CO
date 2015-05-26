import inspect
import time

import numpy as np
from scipy.ndimage.interpolation import map_coordinates
from astropy import units as u
from astropy import log
from scipy import stats

from h2co_modeling import grid_fitter

class generic_paraH2COmodel(object):
    def grid_getmatch_321to303(self, ratio, eratio):
            match,indbest,chi2r = grid_fitter.grid_getmatch(ratio, eratio,
                                                            self.modelratio1)
            return chi2r

    def grid_getmatch_322to321(self, ratio, eratio):
            match,indbest,chi2r = grid_fitter.grid_getmatch(ratio, eratio,
                                                            self.modelratio2)
            return chi2r

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

        h2fromh2co = self.columnarr + np.log10(np.sqrt(np.pi) * linewidth) - h2coabundance
        chi2_h2 = ((h2fromh2co-logh2column)/elogh2column)**2

        return chi2_h2

    def chi2_abundance(self, logabundance, elogabundance):
        # Should there be a factor of np.log10(np.sqrt(np.pi) * linewidth) here?
        # Should linewidth be treated as FWHM instead of sigma, as it presently is?
        model_logabundance = self.columnarr - np.log10(u.pc.to(u.cm)) - self.densityarr
        chi2X = ((model_logabundance-logabundance)/elogabundance)**2
        return chi2X


    def get_parconstraints(self,
                           chi2level=stats.chi2.ppf(stats.norm.cdf(1)-stats.norm.cdf(-1), 3)):
        """
        If parameter constraints have been set with set_constraints or
        set_constraints_fromrow

        Parameters
        ----------
        chi2level : float
            The maximum Delta-chi^2 value to include: this will be used to
            determine the min/max 1-sigma errorbars
        """
        if not hasattr(self, 'chi2'):
            raise AttributeError("Run set_constraints first")

        row = {}

        indbest = np.argmin(self.chi2)
        deltachi2b = (self.chi2-self.chi2.min())
        for parname,pararr in zip(('temperature','column','density'),
                                  (self.temparr,self.columnarr,self.densityarr)):
            row['{0}_chi2'.format(parname)] = pararr.flat[indbest]
            OK = deltachi2b<chi2level
            if np.count_nonzero(OK) > 0:
                row['{0:1.1s}min1sig_chi2'.format(parname)] = pararr[OK].min()
                row['{0:1.1s}max1sig_chi2'.format(parname)] = pararr[OK].max()
            else:
                row['{0:1.1s}min1sig_chi2'.format(parname)] = np.nan
                row['{0:1.1s}max1sig_chi2'.format(parname)] = np.nan

        for parname in ('logh2column', 'elogh2column', 'logabundance',
                        'elogabundance'):
            row[parname] = getattr(self, parname)

        self._parconstraints = row

        return row



    def denstemplot(self):
        self.parplot('dens','tem')

    def denscolplot(self):
        self.parplot('col','dens')

    def coltemplot(self):
        self.parplot('col','tem')
