import numpy as np
import paths
import pyregion
from astropy import coordinates
from astropy import units as u
from astropy import table
from astropy.table import Table,Column
import latex_info
from latex_info import latexdict, exp_to_tex, format_float

tbl = Table.read(paths.tpath('fitted_line_parameters_Chi2Constraints.ipac'), format='ascii.ipac')

def make_column(colname, errcolname, outcolname, unit, formatstr="${0:0.2f}\pm{1:0.2f}$"):
    data = [formatstr.format(d,e) for d,e in zip(tbl[colname], tbl[errcolname])]
    return Column(data=data, name=outcolname, unit=unit)

def make_column_asymm(colname, errlowcolname, errhighcolname, outcolname, unit, formatstr="${0:0.2f}^{{+{1:0.2f}}}_{{-{2:0.2f}}}$"):
    data = [formatstr.format(d,el,eh) for d,el,eh in zip(tbl[colname],
                                                         tbl[colname]-tbl[errlowcolname],
                                                         tbl[errhighcolname]-tbl[colname])]
    return Column(data=data, name=outcolname, unit=unit)

columns = {'Source_Name': 'Source Name',
           'tkin_turb': '$T_{gas, turb}$',
          }
           #'spline_h2coratio321303': '$R_1$',
           #'espline_h2coratio321303': '$\sigma(R_1)$',}

formats = {'tkin_turb': format_float,
          }

outtbl = Table([tbl['Source_Name'],
                make_column('spline_h2coratio321303', 'espline_h2coratio321303', '$R_1$', None),
                make_column('spline_width', 'espline_width', '$\sigma_v$', u.km/u.s),
                make_column('spline_center', 'espline_center', '$v_{lsr}$', u.km/u.s),
                make_column('spline_ampH2CO', 'espline_ampH2CO', '$T_B(H_2CO)$', u.K),
                make_column('logh2column', 'elogh2column', 'log($n(H_2)$)', u.cm**-2),
                make_column_asymm('temperature_chi2', 'tmin1sig_chi2', 'tmax1sig_chi2', '$T_{gas}$', u.K),
                tbl['tkin_turb'],
               ]
              )



latexdict['header_start'] = '\label{tab:observations}'
latexdict['caption'] = '\\formaldehyde Line Parameters and Fit Properties'
latexdict['tablefoot'] = ('\par\n')
#latexdict['col_align'] = 'lllrr'
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
outtbl.write(paths.tpath('line_props.tex'), format='ascii.latex',
             latexdict=latexdict,
             formats=formats,
            )
 

