import numpy as np
import paths
import pyregion
from astropy import coordinates
from astropy import units as u
from astropy import table
from astropy.table import Table,Column
import latex_info
from latex_info import latexdict, exp_to_tex, format_float

tbl = Table.read(paths.tpath('PPV_H2CO_Temperature.ipac'), format='ascii.ipac')

def make_column(colname, errcolname, outcolname, unit, formatstr="${0:0.2f}\pm{1:0.2f}$"):
    data = [formatstr.format(d,e) for d,e in zip(tbl[colname], tbl[errcolname])]
    return Column(data=data, name=outcolname, unit=unit)

def make_column_asymm(colname, errlowcolname, errhighcolname, outcolname, unit, formatstr="${0:0.2f}^{{+{1:0.2f}}}_{{-{2:0.2f}}}$"):
    data = [formatstr.format(d,el,eh) for d,el,eh in zip(tbl[colname],
                                                         tbl[colname]-tbl[errlowcolname],
                                                         tbl[errhighcolname]-tbl[colname])]
    return Column(data=data, name=outcolname, unit=unit)

columns = {'_idx': 'Source ID',
           'DespoticTem': '$T_{gas, turb}$',
           'logh2column': 'log($n(H_2)$)',
          }
           #'spline_h2coratio321303': '$R_1$',
           #'espline_h2coratio321303': '$\sigma(R_1)$',}

def format_float(st):
    return exp_to_tex("{0:0.3g}".format(st))
def format_int(st):
    return ("{0:d}".format(int(st)))

formats = {'DespoticTem': format_int,
           'logh2column': format_float,
           'v_rms': format_float,
           'v_cen': format_float,
           '$\sigma_v$': format_float,
           '$v_{lsr}$':  format_float,
           'Max $T_B(3_{0,3})$': format_float,
           #'$T_{gas}$': format_float,
          }

outtbl = Table([tbl['_idx'],
                make_column('ratio321303', 'eratio321303', '$R_1$', None, "${0:0.3f}\pm{1:0.3f}$"),
                Column(tbl['v_rms']/1e3, name='$\sigma_v$', unit=u.km/u.s),
                Column(tbl['v_cen']/1e3, name='$v_{lsr}$', unit=u.km/u.s),
                Column(tbl['Smax303'], name='Max $T_B(3_{0,3})$', unit=u.K),
                #make_column('spline_ampH2CO', 'espline_ampH2CO', '$T_B(H_2CO)$', u.K),
                #make_column('logh2column', 'elogh2column', 'log($n(H_2)$)', u.cm**-2),
                tbl['logh2column'],
                make_column_asymm('temperature_chi2', 'tmin1sig_chi2', 'tmax1sig_chi2', '$T_{gas}$', u.K),
                Column(data=tbl['DespoticTem'], name='DespoticTem', unit=u.K),
               ]
              )

for old, new in columns.items():
    outtbl.rename_column(old, new)
    if old in formats:
        formats[new] = formats[old]


latexdict['header_start'] = '\label{tab:dendroregions}'
latexdict['caption'] = '\\formaldehyde Parameters and Fit Properties for dendrogram-selected clumps'
latexdict['tablefoot'] = ('\par\n')
#latexdict['col_align'] = 'lllrr'
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
outtbl.write(paths.tpath('dendro_props.tex'), format='ascii.latex',
             latexdict=latexdict,
             formats=formats,
            )
 

outtbl[::10].write(paths.tpath('dendro_props_excerpt.tex'), format='ascii.latex',
             latexdict=latexdict,
             formats=formats,
            )
