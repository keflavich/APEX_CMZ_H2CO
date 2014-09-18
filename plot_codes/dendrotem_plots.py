import numpy as np
import pylab as pl

from paths import hpath
from astropy.table import Table

catalog = Table.read(hpath('PPV_H2CO_Temperature.ipac'), format='ascii.ipac')

sn = (catalog['ratio303321']/catalog['eratio303321'])
sngt50 = sn > 50
ok = np.isfinite(sn)

for ii in range(1,7): pl.figure(ii).clf()

fig1, ax1 = pl.subplots(num=1)
for mask,color in zip((sngt50 & ok, ok & ~sngt50),('b','r')):
    ax1.errorbar(catalog['area_exact'][mask], catalog['temperature_chi2'][mask],
                yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=0.5, marker='.', color=color)
    ax1.set_xscale('log')

fig2, ax2 = pl.subplots(num=2)
for mask,color in zip((ok & sngt50, ok & ~sngt50),('b','r')):
    ax2.errorbar(catalog['ratio303321'][mask], catalog['temperature_chi2'][mask],
                 yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                 #xerr=[catalog['eratio303321'][mask], catalog['eratio303321'][mask]],
                 linestyle='none', capsize=0, alpha=0.2, marker='.', color=color, linewidth=0.3)

fig3, ax3 = pl.subplots(num=3)
ax3.hist(sn[sn==sn], bins=50)

fig4, ax4 = pl.subplots(num=4)
for mask,color in zip((ok & sngt50, ok & ~sngt50),('b','r')):
    ax4.errorbar(catalog['density_chi2'][mask], catalog['temperature_chi2'][mask],
                 yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                 xerr=[catalog['elo_d'][mask], catalog['ehi_d'][mask]],
                 linestyle='none', capsize=0, alpha=0.2, marker='.', color=color, linewidth=0.1)

fig5, ax5 = pl.subplots(num=5)
for mask,color in zip((sngt50 & ok, ok & ~sngt50),('b','r')):
    ax5.errorbar(catalog['x_cen'][mask], catalog['temperature_chi2'][mask],
                yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=0.5, marker='.', color=color)

fig6, ax6 = pl.subplots(num=6)
for mask,color in zip((sngt50 & ok, ok & ~sngt50),('b','r')):
    ax6.errorbar(catalog['v_cen'][mask], catalog['temperature_chi2'][mask],
                yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=0.5, marker='.', color=color)

pl.show()
