import inspect
import time

import numpy as np
from scipy.ndimage.interpolation import map_coordinates
from astropy import units as u
from astropy import log
from scipy import stats

from h2co_modeling import grid_fitter

class generic_paraH2COmodel(object):
    def grid_getmatch_321to303(self, ratio, eratio, chi2_thresh=1):
            match,indbest,chi2r = grid_fitter.grid_getmatch(ratio, eratio,
                                                            self.modelratio1,
                                                            chi2_thresh=chi2_thresh)
            return chi2r

    def grid_getmatch_322to321(self, ratio, eratio, chi2_thresh=1):
            match,indbest,chi2r = grid_fitter.grid_getmatch(ratio, eratio,
                                                            self.modelratio2,
                                                            chi2_thresh=chi2_thresh)
            return chi2r

    @property
    def chi2(self):
        return self._chi2

    @chi2.setter
    def chi2(self, value):
        self._chi2 = value
        self._likelihood = np.exp(-value/2)

    @property
    def likelihood(self):
        return self._likelihood

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
                           nsigma=1):
        """
        If parameter constraints have been set with set_constraints or
        set_constraints_fromrow

        Parameters
        ----------
        nsigma : float
            The number of sigmas to go out to when determining errors
        """
        if not hasattr(self, 'chi2'):
            raise AttributeError("Run set_constraints first")

        row = {}

        inds = np.argsort(self.likelihood.flat)
        cdf = np.cumsum(self.likelihood.flat[inds])
        cdf = cdf/cdf[-1]
        frac_above = (stats.norm.cdf(nsigma)-stats.norm.cdf(-nsigma))
        cdfmin = np.argmin(np.abs(cdf - (1-frac_above)))
        sigma_like = self.likelihood.flat[inds][cdfmin]

        indbest = np.argmax(self.likelihood)
        for parname,pararr in zip(('temperature','column','density'),
                                  (self.temparr,self.columnarr,self.densityarr)):
            row['{0}_chi2'.format(parname)] = pararr.flat[indbest]
            row['expected_{0}'.format(parname)] = (pararr*self.likelihood).sum() / self.likelihood.sum()
            OK = self.likelihood > sigma_like
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

    @property
    def parconstraints(self):
        if not hasattr(self,'_parconstraints') or self._parconstraints is None:
            return self.get_parconstraints()
        else:
            return self._parconstraints


    def denstemplot(self):
        self.parplot('dens','tem')

    def denscolplot(self):
        self.parplot('col','dens')

    def coltemplot(self):
        self.parplot('col','tem')
