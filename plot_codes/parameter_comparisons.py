import os
import pylab as pl
import paths
from astropy import table
from paths import analysispath
import numpy as np
from astropy import coordinates
from astropy import units as u
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

pcfittable = table.Table.read(os.path.join(analysispath,
                                         'fitted_line_parameters_Chi2Constraints.ipac'),
                            format='ascii.ipac')

lolim = pcfittable['tmax1sig_chi2'] > 340
maps = np.char.startswith(pcfittable['Source_Name'], 'Map')
ok = ~np.isnan(pcfittable['tmin1sig_chi2']) & (pcfittable['width'] < 40)
flags = {'is_map': maps,
         'is_lolim': lolim,
         'is_ok': ok}
# Don't plot these for now...
pcfittable = pcfittable[(~lolim) & ok]
maps = np.char.startswith(pcfittable['Source_Name'], 'Map')

fig4 = pl.figure(4)
fig4.clf()
ax = fig4.add_subplot(1,3,1)
ax.errorbar(pcfittable['temperature_chi2'], pcfittable['density_chi2'],
            yerr=[pcfittable['density_chi2']-pcfittable['dmin1sig_chi2'],
                  pcfittable['dmax1sig_chi2']-pcfittable['density_chi2']],
            xerr=[pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'],
                  pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2']],
            linestyle='none', marker='s', linewidth=1, alpha=0.5)
ax2 = fig4.add_subplot(1,3,2)
# Don't do this any more: it relies on having the RADEX fits, which we don't.
#ax2.errorbar(pcfittable['temperature_chi2'], pcfittable['temperature'],
#             yerr=pcfittable['etemperature'],
#             xerr=[pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'],
#                   pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2']],
#             linestyle='none', marker='s', linewidth=1, alpha=0.5)
ax2.plot([0,300],[0,300],'k--',linewidth=2,alpha=0.5)

fig5 = pl.figure(5)
fig5.clf()
ax5 = fig5.gca()
ax5.errorbar(coordinates.Angle(pcfittable['GLON']*u.deg).wrap_at(180*u.deg).value[maps],
             pcfittable['temperature_chi2'][maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax5.set_ylim(0,150)
ax5.set_ylabel("Kinetic Temperature (K)")
ax5.set_xlabel("Galactic Longitude ($^{\\circ}$)")
fig5.savefig(paths.fpath('chi2_temperature_vs_glon_byfield.pdf'),
                         bbox_inches='tight')
ax5.errorbar(coordinates.Angle(pcfittable['GLON']*u.deg).wrap_at(180*u.deg).value[~maps],
             pcfittable['temperature_chi2'][~maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[~maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5)
fig5.savefig(paths.fpath('chi2_temperature_vs_glon_fieldsandsources.pdf'),
                         bbox_inches='tight')

fig6 = pl.figure(6)
fig6.clf()
ax6 = fig6.gca()
ax6.errorbar(pcfittable['higaldusttem'][maps],
             pcfittable['temperature_chi2'][maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax6.plot([15,30],[15,30],'k--')
ax6.set_xlabel("HiGal Fitted Temperature")
ax6.set_ylabel("Kinetic Temperature (K)")
ax6.set_ylim(0,200)
ax6.set_xlim(15,30)
fig6.savefig(paths.fpath('chi2_temperature_vs_higaltemperature_byfield.pdf'),
                         bbox_inches='tight')
ax6.errorbar(pcfittable['higaldusttem'][~maps],
             pcfittable['temperature_chi2'][~maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[~maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
ax6.set_ylim(0,200)
ax6.set_xlim(15,30)
fig6.savefig(paths.fpath('chi2_temperature_vs_higaltemperature_fieldsandsources.pdf'),
                         bbox_inches='tight')

fig7 = pl.figure(7)
fig7.clf()
ax7 = fig7.gca()
ax7.errorbar(pcfittable['width'][maps]*(8*np.log(2))**0.5,
             pcfittable['temperature_chi2'][maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax7.set_xlabel("Line FWHM (km s$^{-1}$)")
ax7.set_ylabel("Kinetic Temperature (K)")
fig7.savefig(paths.fpath('chi2_temperature_vs_linewidth_byfield.pdf'),
                         bbox_inches='tight')
ax7.errorbar(pcfittable['width'][~maps]*(8*np.log(2))**0.5,
             pcfittable['temperature_chi2'][~maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[~maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig7.savefig(paths.fpath('chi2_temperature_vs_linewidth_fieldsandsources.pdf'),
                         bbox_inches='tight')


fig8 = pl.figure(8)
fig8.clf()
ax8 = fig8.gca()
ax8.errorbar(pcfittable['ampH2CO'][maps],
             pcfittable['temperature_chi2'][maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax8.set_xlabel("H2CO Peak Amplitude")
ax8.set_ylabel("Kinetic Temperature (K)")
fig8.savefig(paths.fpath('chi2_temperature_vs_h2coamp_byfield.pdf'),
                         bbox_inches='tight')
ax8.errorbar(pcfittable['ampH2CO'][~maps],
             pcfittable['temperature_chi2'][~maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[~maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig8.savefig(paths.fpath('chi2_temperature_vs_h2coamp_fieldsandsources.pdf'),
                         bbox_inches='tight')


fig9 = pl.figure(9)
fig9.clf()
ax9 = fig9.gca()
ax9.set_xscale('log')
ax9.errorbar(pcfittable['higalcolumndens'][maps],
             pcfittable['temperature_chi2'][maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax9.set_xlabel("Hi-Gal Fitted Column Density")
ax9.set_ylabel("Kinetic Temperature (K)")
fig9.savefig(paths.fpath('chi2_temperature_vs_higalcolumn_byfield.pdf'),
                         bbox_inches='tight')
ax9.errorbar(pcfittable['higalcolumndens'][~maps],
             pcfittable['temperature_chi2'][~maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[~maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig9.savefig(paths.fpath('chi2_temperature_vs_higalcolumn_fieldsandsources.pdf'),
                         bbox_inches='tight')

fig10 = pl.figure(10)
fig10.clf()
ax10 = fig10.gca()
ax10.errorbar(pcfittable['width'][maps]*(8*np.log(2))**0.5,
             pcfittable['h2coratio321303'][maps],
             yerr=pcfittable['eh2coratio321303'][maps],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax10.set_xlabel("Line FWHM (km s$^{-1}$)")
ax10.set_ylabel("Ratio 321/303")
fig10.savefig(paths.fpath('ratio_vs_linewidth_byfield.pdf'),
                         bbox_inches='tight')
ax10.errorbar(pcfittable['width'][~maps]*(8*np.log(2))**0.5,
              pcfittable['h2coratio321303'][~maps],
              yerr=pcfittable['eh2coratio321303'][~maps],
              linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig10.savefig(paths.fpath('ratio_vs_linewidth_fieldsandsources.pdf'),
                         bbox_inches='tight')


fig11 = pl.figure(11)
fig11.clf()
ax11 = fig11.gca()
ax11.errorbar(pcfittable['higaldusttem'][maps],
             pcfittable['h2coratio321303'][maps],
             yerr=pcfittable['eh2coratio321303'][maps],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax11.set_ylim(0,200)
ax11.set_xlim(15,30)
ax11.set_xlabel("HiGal Fitted Temperature")
ax11.set_ylabel("Ratio 321/303")
fig11.savefig(paths.fpath('ratio_vs_higaltemperature_byfield.pdf'),
                         bbox_inches='tight')
ax11.errorbar(pcfittable['higaldusttem'][~maps],
              pcfittable['h2coratio321303'][~maps],
              yerr=pcfittable['eh2coratio321303'][~maps],
              linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig11.savefig(paths.fpath('ratio_vs_higaltemperature_fieldsandsources.pdf'),
                         bbox_inches='tight')

# RADEX fitting has been removed
#fig12 = pl.figure(12)
#fig12.clf()
#ax = fig12.add_subplot(1,1,1)
#ax.errorbar(pcfittable['temperature_chi2'], pcfittable['temperature'],
#            yerr=pcfittable['etemperature'],
#            xerr=[pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'],
#                  pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2']],
#            linestyle='none', marker='s', linewidth=1, alpha=0.5)
#ax.plot([0,300],[0,300],'k--',linewidth=2,alpha=0.5)
#ax.set_title("DEBUG: RADEX+pyspeckit-fitted temperature vs. $\\chi^2$ temperature")
#ax.set_xlabel("$\\chi^2$ Temperature")
#ax.set_ylabel("RADEX+pyspeckit Temperature")
#ax.axis([0,350,0,350])


fig13 = pl.figure(13)
fig13.clf()
ax13 = fig13.gca()
ax13.errorbar(pcfittable['area'][maps],
              pcfittable['temperature_chi2'][maps],
              yerr=[pcfittable['temperature_chi2'][maps]-pcfittable['tmin1sig_chi2'][maps],
                    pcfittable['tmax1sig_chi2'][maps]-pcfittable['temperature_chi2'][maps]],
              linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax13.set_xlabel("Area (square degrees)")
ax13.set_ylabel("Kinetic Temperature (K)")
ax13.set_xscale('log')
fig13.savefig(paths.fpath('temperature_vs_area_byfield.pdf'),
                         bbox_inches='tight')
ax13.errorbar(pcfittable['area'][~maps],
              pcfittable['temperature_chi2'][~maps],
              yerr=[pcfittable['temperature_chi2'][~maps]-pcfittable['tmin1sig_chi2'][~maps],
                    pcfittable['tmax1sig_chi2'][~maps]-pcfittable['temperature_chi2'][~maps]],
              linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig13.savefig(paths.fpath('temperature_vs_area_fieldsandsources.pdf'),
                         bbox_inches='tight')


fig14 = pl.figure(14)
fig14.clf()
ax14 = fig14.gca()
ax14.errorbar(pcfittable['higalcolumndens'][maps],
             pcfittable['temperature_chi2'][maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
#ax14.plot([15,30],[15,30],'k--')
ax14.set_xlabel("HiGal Fitted Column Density")
ax14.set_ylabel("Kinetic Temperature (K)")
fig14.savefig(paths.fpath('chi2_temperature_vs_higaldustcol_byfield.pdf'),
                         bbox_inches='tight')
ax14.errorbar(pcfittable['higalcolumndens'][~maps],
             pcfittable['temperature_chi2'][~maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[~maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig14.savefig(paths.fpath('chi2_temperature_vs_higaldustcol_fieldsandsources.pdf'),
                         bbox_inches='tight')

# pcfittable[np.abs(pcfittable['temperature_chi2']-pcfittable['higaldusttem'])/pcfittable['higaldusttem'] < 1.5].pprint()
