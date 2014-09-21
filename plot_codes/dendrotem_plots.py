import numpy as np
import pylab as pl

from paths import hpath
from astropy.table import Table

catalog = Table.read(hpath('PPV_H2CO_Temperature.ipac'), format='ascii.ipac')

sn = (catalog['ratio303321']/catalog['eratio303321'])
sngt50 = sn > 50
sn25_50 = (sn > 25) & (sn < 50)
ok = np.isfinite(sn) & (sn>5)

for ii in range(1,15): pl.figure(ii).clf()

masks_colors = zip((ok & ~sngt50 & ~sn25_50, sn25_50 & ok, sngt50 & ok, ),
                   ('b','g','r'),
                   (0.2,0.3,0.4),
                  )

fig1, ax1 = pl.subplots(num=1)
for mask,color,alpha in masks_colors:
    ax1.errorbar(catalog['area_exact'][mask], catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax1.set_xscale('log')
    ax1.set_xlabel("Size Scale (pixels?)")
    ax1.set_ylabel("Temperature")

fig2, ax2 = pl.subplots(num=2)
for mask,color,alpha in masks_colors:
    ax2.errorbar(catalog['ratio303321'][mask], catalog['temperature_chi2'][mask],
                 #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                 #xerr=[catalog['eratio303321'][mask], catalog['eratio303321'][mask]],
                 linestyle='none', capsize=0, alpha=alpha, marker='.', color=color, linewidth=0.3)

fig3, ax3 = pl.subplots(num=3)
ax3.hist(sn[sn==sn], bins=50)

fig4, ax4 = pl.subplots(num=4)
for mask,color,alpha in masks_colors:
    ax4.errorbar(catalog['density_chi2'][mask], catalog['temperature_chi2'][mask],
                 #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                 #xerr=[catalog['elo_d'][mask], catalog['ehi_d'][mask]],
                 linestyle='none', capsize=0, alpha=alpha, marker='.', color=color, linewidth=0.1)

fig5, ax5 = pl.subplots(num=5)
for mask,color,alpha in masks_colors:
    lon = catalog['x_cen'][mask]
    lon[lon>180] -= 360
    ax5.errorbar(lon, catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)

fig6, ax6 = pl.subplots(num=6)
for mask,color,alpha in masks_colors:
    ax6.errorbar(catalog['v_cen'][mask], catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax6.set_xlabel("Centroid Velocity")

fig7, ax7 = pl.subplots(num=7)
for mask,color,alpha in masks_colors:
    ax7.errorbar(catalog['radius'][mask], catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax7.set_xscale('log')
    ax7.set_xlabel("Radius (pixels?)")
    ax7.set_ylabel("Temperature")

fig8, ax8 = pl.subplots(num=8)
for mask,color,alpha in masks_colors:
    ax8.errorbar(catalog['Smean303'][mask], catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                #xerr=[catalog['e303'][mask], catalog['e303'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax8.set_xlabel("$S(3_{0,3}-2_{0,2})$ (K)")
    ax8.set_ylabel("Temperature (K)")

fig9, ax9 = pl.subplots(num=9)
for mask,color,alpha in masks_colors:
    ax9.errorbar(catalog['Smean321'][mask], catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                #xerr=[catalog['e321'][mask], catalog['e321'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax9.set_xlabel("$S(3_{2,1}-2_{2,0})$ (K)")
    ax9.set_ylabel("Temperature (K)")

fig10, ax10 = pl.subplots(num=10)
for mask,color,alpha in masks_colors:
    ax10.errorbar(catalog['13comean'][mask], catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax10.set_xlabel("${S}(^{13}$CO) (K)")
    ax10.set_ylabel("Temperature (K)")

fig11, ax11 = pl.subplots(num=11)
for mask,color,alpha in masks_colors:
    ax11.errorbar(catalog['13comean'][mask], catalog['Smean303'][mask],
                #yerr=[catalog['e303'][mask], catalog['e303'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax11.set_xlabel("$\\bar{S}(^{13}$CO) (K)")
    ax11.set_ylabel("$S(3_{0,3}-2_{0,2})$ (K)")

fig12, ax12 = pl.subplots(num=12)
for mask,color,alpha in masks_colors:
    ax12.errorbar(catalog['v_rms'][mask], catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax12.set_xlabel("RMS Velocity")

pl.show()
