import matplotlib
import paths
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))
import os
import pylab as pl
from astropy import table
from paths import analysispath
import numpy as np
from astropy import coordinates
from astropy import units as u
import heating

pcfittable = table.Table.read(os.path.join(analysispath,
                                         'fitted_line_parameters_Chi2Constraints.ipac'),
                            format='ascii.ipac')

lolim = pcfittable['tmax1sig_chi2'] > 340
maps = np.char.startswith(pcfittable['Source_Name'], 'Map')
ok = ~np.isnan(pcfittable['tmin1sig_chi2']) & (pcfittable['width'] < 40) & (pcfittable['h2coratio321303']/pcfittable['eh2coratio321303'] > 5) & pcfittable['is_good'].astype('bool')
flags = {'is_map': maps,
         'is_lolim': lolim,
         'is_ok': ok}
# Don't plot these for now...
pcfittable = pcfittable[(~lolim) & ok]
maps = np.char.startswith(pcfittable['Source_Name'], 'Map')
lolim_conservative = pcfittable['tmax1sig_chi2'] > 150

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
             capsize=0, markeredgecolor='none',
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax5.set_ylim(0,150)
ax5.set_ylabel("Temperature (K)")
ax5.set_xlabel("Galactic Longitude ($^{\\circ}$)")
fig5.savefig(paths.fpath('chi2_temperature_vs_glon_byfield.pdf'),
                         bbox_inches='tight')
ax5.errorbar(coordinates.Angle(pcfittable['GLON']*u.deg).wrap_at(180*u.deg).value[~maps],
             pcfittable['temperature_chi2'][~maps],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[~maps],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[~maps]],
             capsize=0, markeredgecolor='none',
             linestyle='none', marker='s', linewidth=1, alpha=0.5)
fig5.savefig(paths.fpath('chi2_temperature_vs_glon_fieldsandsources.pdf'),
                         bbox_inches='tight')

fig6 = pl.figure(6)
fig6.clf()
ax6 = fig6.gca()
mask = maps&~lolim_conservative
ax6.errorbar(pcfittable['higaldusttem'][mask],
             pcfittable['temperature_chi2'][mask],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[mask],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[mask]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r', capsize=0)
ax6.plot([15,30],[15,30],'k--')
mask = maps&lolim_conservative
ax6.plot(pcfittable['higaldusttem'][mask],
         pcfittable['tmin1sig_chi2'][mask],
         marker='^',
         markersize=10,
         markeredgecolor='none',
         color='r',
         alpha=0.5,
         linestyle='none')
ax6.set_xlabel("HiGal Dust Temperature (K)")
ax6.set_ylabel("H$_2$CO Temperature (K)")
ax6.set_ylim(0,200)
ax6.set_xlim(15,30)
fig6.savefig(paths.fpath('chi2_temperature_vs_higaltemperature_byfield.pdf'),
                         bbox_inches='tight')
mask = (~maps)&(~lolim_conservative)
ax6.errorbar(pcfittable['higaldusttem'][mask],
             pcfittable['temperature_chi2'][mask],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[mask],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[mask]],
             capsize=0,
             markeredgecolor='none',
             markersize=10,
             linestyle='none', marker='s', linewidth=0.5, alpha=0.5, color='b')

mask = (~maps)&lolim_conservative
ax6.plot(pcfittable['higaldusttem'][mask],
         pcfittable['tmin1sig_chi2'][mask],
         marker='^',
         markersize=10,
         markeredgecolor='none',
         color='b',
         alpha=0.5,
         linestyle='none')

ax6.set_ylim(10,150)
ax6.set_xlim(15,30)
fig6.savefig(paths.fpath('chi2_temperature_vs_higaltemperature_fieldsandsources_notitle.pdf'),
                         bbox_inches='tight')
ax6.set_title("Hand-selected regions")
fig6.savefig(paths.fpath('chi2_temperature_vs_higaltemperature_fieldsandsources.pdf'),
                         bbox_inches='tight')

fig7 = pl.figure(7)
fig7.clf()
ax7 = fig7.gca()
mask = maps&~lolim_conservative
ax7.errorbar(pcfittable['width'][mask]*(8*np.log(2))**0.5,
             pcfittable['temperature_chi2'][mask],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[mask],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[mask]],
             capsize=0,
             markersize=10,
             markeredgecolor='none',
             linestyle='none', marker='s', linewidth=0.5, alpha=0.6, color='r')
mask = maps&lolim_conservative
ax7.plot(pcfittable['width'][mask]*(8*np.log(2))**0.5,
         pcfittable['tmin1sig_chi2'][mask],
         marker='^',
         markersize=10,
         markeredgecolor='none',
         color='r',
         alpha=0.4,
         linestyle='none')

linewidths = np.linspace(0,pcfittable['width'].max())*u.km/u.s
ax7.plot(linewidths*2.35, [heating.tkin_all(10**4*u.cm**-3, sigma, 10*u.pc,
                                            5*u.km/u.s/u.pc, 30*u.K)
                           for sigma in linewidths],
        linestyle='--', color='k', label='$n=10^4$ cm$^{-3}$', zorder=-5)
ax7.plot(linewidths*2.35, [heating.tkin_all(10**4*u.cm**-3, sigma, 10*u.pc,
                                            1*u.km/u.s/u.pc, 30*u.K)
                           for sigma in linewidths],
        linestyle='--', color='r', label='$n=10^4$ cm$^{-3}$, $dv/dr=1$', zorder=-5, linewidth=2, alpha=0.5)
ax7.plot(linewidths*2.35, [heating.tkin_all(10**4*u.cm**-3, sigma, 20*u.pc,
                                            5*u.km/u.s/u.pc, 30*u.K)
                           for sigma in linewidths],
         linestyle='--', color='b', label='$n=10^4$ cm$^{-3}$, $L=20$ pc', zorder=-5, alpha=0.5, linewidth=2)

ax7.plot(linewidths*2.35, [heating.tkin_all(10**5*u.cm**-3, sigma, 10*u.pc,
                                            5*u.km/u.s/u.pc, 30*u.K)
                           for sigma in linewidths],
         linestyle=':', color='k', label='$n=10^5$ cm$^{-3}$', zorder=-5)
ax7.plot(linewidths*2.35, [heating.tkin_all(10**6*u.cm**-3, sigma, 10*u.pc,
                                            5*u.km/u.s/u.pc, 30*u.K)
                           for sigma in linewidths],
        linestyle='-.', color='k', label='$n=10^6$ cm$^{-3}$', zorder=-5)
ax7.plot(linewidths*2.35, [heating.tkin_all(10**5*u.cm**-3, sigma, 10*u.pc,
                                            5*u.km/u.s/u.pc, 30*u.K, crir=1e-15*u.s**-1)
                           for sigma in linewidths],
         linestyle='-', color='g', label='$n=10^5$ cm$^{-3}$, $\zeta_{CR}=10^{-15}$ s$^{-1}$', zorder=-10, alpha=0.25, linewidth=4)
ax7.plot(linewidths*2.35, [heating.tkin_all(10**5*u.cm**-3, sigma, 10*u.pc,
                                            5*u.km/u.s/u.pc, 30*u.K, crir=1e-14*u.s**-1)
                           for sigma in linewidths],
         linestyle=':', color='purple', label='$n=10^5$ cm$^{-3}$, $\zeta_{CR}=10^{-14}$ s$^{-1}$', zorder=-10, alpha=0.25, linewidth=4)

box = ax7.get_position()
ax7.set_position([box.x0, box.y0, box.width * 0.7, box.height])
ax7.legend(loc='center left', fontsize=16, bbox_to_anchor=(1.0, 0.75))
ax7.set_xlabel("Line FWHM (km s$^{-1}$)")
ax7.set_ylabel("Temperature (K)")
ax7.set_ylim(10,150)
fig7.savefig(paths.fpath('chi2_temperature_vs_linewidth_byfield.pdf'),
                         bbox_inches='tight')

mask = (~maps)&(~lolim_conservative)
ax7.errorbar(pcfittable['width'][mask]*(8*np.log(2))**0.5,
             pcfittable['temperature_chi2'][mask],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[mask],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[mask]],
             capsize=0,
             markeredgecolor='none',
             markersize=10,
             linestyle='none', marker='s', linewidth=0.5, alpha=0.6, color='b')

mask = (~maps)&lolim_conservative
ax7.plot(pcfittable['width'][mask]*(8*np.log(2))**0.5,
         pcfittable['tmin1sig_chi2'][mask],
         marker='^',
         markersize=10,
         markeredgecolor='none',
         color='b',
         alpha=0.4,
         linestyle='none')

ax7.set_ylim(10,150)
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
ax8.set_ylabel("Temperature (K)")
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
ax9.set_ylabel("Temperature (K)")
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
ax13.set_ylabel("Temperature (K)")
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
ax14.set_ylabel("Temperature (K)")
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

fig15 = pl.figure(15)
fig15.clf()
ax15 = fig15.gca()
mask = maps&~lolim_conservative
ax15.errorbar(pcfittable['tkin_turb'][mask],
              pcfittable['temperature_chi2'][mask],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[mask],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[mask]],
             capsize=0,
             markersize=10,
             markeredgecolor='none',
             linestyle='none', marker='s', linewidth=0.5, alpha=0.6, color='r')
mask = maps&lolim_conservative
ax15.plot(pcfittable['tkin_turb'][mask],
         pcfittable['tmin1sig_chi2'][mask],
         marker='^',
         markersize=10,
         markeredgecolor='none',
         color='r',
         alpha=0.4,
         linestyle='none')

mask = (maps) & (~lolim_conservative) & ((pcfittable['tmin1sig_chi2'] > pcfittable['tkin_turb']) | (pcfittable['tmax1sig_chi2'] < pcfittable['tkin_turb']))
ax15.plot(pcfittable['tkin_turb'][mask],
         pcfittable['temperature_chi2'][mask],
         marker='s',
         markersize=15,
         markeredgecolor='r',
         markerfacecolor='none',
         markeredgewidth=0.5,
         alpha=0.4,
         linestyle='none')
mask = (maps) & (lolim_conservative) & ((pcfittable['tmin1sig_chi2'] > pcfittable['tkin_turb']))
ax15.plot(pcfittable['tkin_turb'][mask],
         pcfittable['tmin1sig_chi2'][mask],
         marker='^',
         markersize=15,
         markeredgecolor='r',
         markerfacecolor='none',
         markeredgewidth=0.5,
         alpha=0.4,
         linestyle='none')

# Sources with T_predicted >> T_measured
#high_badpredictions = (pcfittable['tkin_turb'] > pcfittable['tmax1sig_chi2'])&(~lolim_conservative)
#high_badpredictions = (pcfittable['tkin_turb'] > 120)&(~lolim_conservative)
#for row,is_map in zip(pcfittable[high_badpredictions], maps[high_badpredictions]):
#    xy = np.array((row['tkin_turb'], row['temperature_chi2']))
#    ax15.annotate("{0}_{1}".format(row['Source_Name'], row['ComponentID']),
#                  xy,
#                  xytext=xy-(15, 7),
#                  color='r' if is_map else 'b'
#                 )

ax15.plot([0,200], [0,200], 'k--', alpha=0.5, zorder=-5)
ax15.set_xlabel("Turbulence-driven Temperature (K)")
ax15.set_ylabel("H$_2$CO Temperature (K)")
ax15.set_ylim(10,150)
ax15.set_xlim(10,180)
fig15.savefig(paths.fpath('chi2_temperature_vs_turbulenttemperature_byfield.pdf'),
                         bbox_inches='tight')

mask = (~maps)&(~lolim_conservative)
ax15.errorbar(pcfittable['tkin_turb'][mask],
             pcfittable['temperature_chi2'][mask],
             yerr=[(pcfittable['temperature_chi2']-pcfittable['tmin1sig_chi2'])[mask],
                   (pcfittable['tmax1sig_chi2']-pcfittable['temperature_chi2'])[mask]],
             capsize=0,
             markeredgecolor='none',
             markersize=10,
             linestyle='none', marker='s', linewidth=0.5, alpha=0.6, color='b')

mask = (~maps)&lolim_conservative
ax15.plot(pcfittable['tkin_turb'][mask],
         pcfittable['tmin1sig_chi2'][mask],
         marker='^',
         markersize=10,
         markeredgecolor='none',
         color='b',
         alpha=0.4,
         linestyle='none')

mask = (~maps) & (~lolim_conservative) & ((pcfittable['tmin1sig_chi2'] > pcfittable['tkin_turb']) | (pcfittable['tmax1sig_chi2'] < pcfittable['tkin_turb']))
ax15.plot(pcfittable['tkin_turb'][mask],
         pcfittable['temperature_chi2'][mask],
         marker='s',
         markersize=15,
         markeredgecolor='b',
         markerfacecolor='none',
         markeredgewidth=0.5,
         alpha=0.4,
         linestyle='none')
mask = (~maps) & (lolim_conservative) & ((pcfittable['tmin1sig_chi2'] > pcfittable['tkin_turb']))
ax15.plot(pcfittable['tkin_turb'][mask],
         pcfittable['tmin1sig_chi2'][mask],
         marker='^',
         markersize=15,
         markeredgecolor='b',
         markerfacecolor='none',
         markeredgewidth=0.5,
         alpha=0.4,
         linestyle='none')

ax15.set_ylim(10,150)
fig15.savefig(paths.fpath('chi2_temperature_vs_turbulenttemperature_fieldsandsources_notitle.pdf'),
                         bbox_inches='tight')
ax15.set_title("Hand-selected regions")
fig15.savefig(paths.fpath('chi2_temperature_vs_turbulenttemperature_fieldsandsources.pdf'),
                         bbox_inches='tight')
