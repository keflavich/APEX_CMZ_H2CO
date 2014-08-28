import os
import pylab as pl
import paths
from astropy import table
from paths import analysispath
import numpy as np
from astropy import coordinates
from astropy import units as u

fittable = table.Table.read(os.path.join(analysispath,
                                         'fitted_line_parameters_Chi2Constraints.ipac'),
                            format='ascii.ipac')

lolim = fittable['tmax1sig_chi2'] > 340
maps = np.char.startswith(fittable['Source_Name'], 'Map')
ok = ~np.isnan(fittable['tmin1sig_chi2']) & (fittable['width'] < 40)
flags = {'is_map': maps,
         'is_lolim': lolim,
         'is_ok': ok}
# Don't plot these for now...
fittable = fittable[(~lolim) & ok]
maps = np.char.startswith(fittable['Source_Name'], 'Map')

fig4 = pl.figure(4)
fig4.clf()
ax = fig4.add_subplot(1,3,1)
ax.errorbar(fittable['temperature_chi2'], fittable['density_chi2'],
            yerr=[fittable['density_chi2']-fittable['dmin1sig_chi2'],
                  fittable['dmax1sig_chi2']-fittable['density_chi2']],
            xerr=[fittable['temperature_chi2']-fittable['tmin1sig_chi2'],
                  fittable['tmax1sig_chi2']-fittable['temperature_chi2']],
            linestyle='none', marker='s', linewidth=1, alpha=0.5)
ax2 = fig4.add_subplot(1,3,2)
ax2.errorbar(fittable['temperature_chi2'], fittable['temperature'],
             yerr=fittable['etemperature'],
             xerr=[fittable['temperature_chi2']-fittable['tmin1sig_chi2'],
                   fittable['tmax1sig_chi2']-fittable['temperature_chi2']],
             linestyle='none', marker='s', linewidth=1, alpha=0.5)
ax2.plot([0,300],[0,300],'k--',linewidth=2,alpha=0.5)

fig5 = pl.figure(5)
fig5.clf()
ax5 = fig5.gca()
ax5.errorbar(coordinates.Angle(fittable['GLON']*u.deg).wrap_at(180*u.deg).value[maps],
             fittable['temperature_chi2'][maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax5.set_ylim(0,150)
ax5.set_ylabel("Kinetic Temperature (K)")
ax5.set_xlabel("Galactic Longitude ($^{\\circ}$)")
fig5.savefig(paths.fpath('chi2_temperature_vs_glon_byfield.pdf'),
                         bbox_inches='tight')
ax5.errorbar(coordinates.Angle(fittable['GLON']*u.deg).wrap_at(180*u.deg).value[~maps],
             fittable['temperature_chi2'][~maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[~maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5)
fig5.savefig(paths.fpath('chi2_temperature_vs_glon_fieldsandsources.pdf'),
                         bbox_inches='tight')

fig6 = pl.figure(6)
fig6.clf()
ax6 = fig6.gca()
ax6.errorbar(fittable['higaldusttem'][maps],
             fittable['temperature_chi2'][maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax6.plot([15,30],[15,30],'k--')
ax6.set_xlabel("HiGal Fitted Temperature")
ax6.set_ylabel("Kinetic Temperature (K)")
fig6.savefig(paths.fpath('chi2_temperature_vs_higaltemperature_byfield.pdf'),
                         bbox_inches='tight')
ax6.errorbar(fittable['higaldusttem'][~maps],
             fittable['temperature_chi2'][~maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[~maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig6.savefig(paths.fpath('chi2_temperature_vs_higaltemperature_fieldsandsources.pdf'),
                         bbox_inches='tight')

fig7 = pl.figure(7)
fig7.clf()
ax7 = fig7.gca()
ax7.errorbar(fittable['width'][maps]*(8*np.log(2))**0.5,
             fittable['temperature_chi2'][maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax7.set_xlabel("Line FWHM (km s$^{-1}$)")
ax7.set_ylabel("Kinetic Temperature (K)")
fig7.savefig(paths.fpath('chi2_temperature_vs_linewidth_byfield.pdf'),
                         bbox_inches='tight')
ax7.errorbar(fittable['width'][~maps]*(8*np.log(2))**0.5,
             fittable['temperature_chi2'][~maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[~maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig7.savefig(paths.fpath('chi2_temperature_vs_linewidth_fieldsandsources.pdf'),
                         bbox_inches='tight')


fig8 = pl.figure(8)
fig8.clf()
ax8 = fig8.gca()
ax8.errorbar(fittable['ampH2CO'][maps],
             fittable['temperature_chi2'][maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax8.set_xlabel("H2CO Peak Amplitude")
ax8.set_ylabel("Kinetic Temperature (K)")
fig8.savefig(paths.fpath('chi2_temperature_vs_h2coamp_byfield.pdf'),
                         bbox_inches='tight')
ax8.errorbar(fittable['ampH2CO'][~maps],
             fittable['temperature_chi2'][~maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[~maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig8.savefig(paths.fpath('chi2_temperature_vs_h2coamp_fieldsandsources.pdf'),
                         bbox_inches='tight')


fig9 = pl.figure(9)
fig9.clf()
ax9 = fig9.gca()
ax9.set_xscale('log')
ax9.errorbar(fittable['higalcolumndens'][maps],
             fittable['temperature_chi2'][maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='r')
ax9.set_xlabel("Hi-Gal Fitted Column Density")
ax9.set_ylabel("Kinetic Temperature (K)")
fig9.savefig(paths.fpath('chi2_temperature_vs_higalcolumn_byfield.pdf'),
                         bbox_inches='tight')
ax9.errorbar(fittable['higalcolumndens'][~maps],
             fittable['temperature_chi2'][~maps],
             yerr=[(fittable['temperature_chi2']-fittable['tmin1sig_chi2'])[~maps],
                   (fittable['tmax1sig_chi2']-fittable['temperature_chi2'])[~maps]],
             linestyle='none', marker='s', linewidth=1, alpha=0.5, color='b')
fig9.savefig(paths.fpath('chi2_temperature_vs_higalcolumn_fieldsandsources.pdf'),
                         bbox_inches='tight')

