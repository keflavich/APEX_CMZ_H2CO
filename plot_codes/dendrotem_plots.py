import numpy as np
import pylab as pl

from paths import hpath, apath
from astropy.table import Table, Column
from astropy import log
from scipy.interpolate import PiecewisePolynomial

catalog = Table.read(hpath('PPV_H2CO_Temperature.ipac'), format='ascii.ipac')

sn = (catalog['ratio303321']/catalog['eratio303321'])
sngt50 = sn > 50
sn25_50 = (sn > 25) & (sn < 50)
ok = np.isfinite(sn) & (sn>5)

for ii in range(1,13): pl.figure(ii).clf()

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
                 linestyle='none', capsize=0, alpha=alpha, marker='.',
                 color=color, linewidth=0.3)
    ax2.set_xlabel("Ratio 303/321")
    ax2.set_ylabel("Temperature")
    
# Determine approximate best-fit
sel = catalog['temperature_chi2'][ok] < 60
fparslt60 = np.polyfit(catalog['ratio303321'][ok][sel],
                       catalog['temperature_chi2'][ok][sel], 2)
p1 = np.polynomial.polynomial.Polynomial(fparslt60[::-1])

sel = ((catalog['temperature_chi2'][ok] > 60) &
       (catalog['temperature_chi2'][ok] < 150))
fparsgt60 = np.polyfit(catalog['ratio303321'][ok][sel],
                       catalog['temperature_chi2'][ok][sel], 2)
p2 = np.polynomial.polynomial.Polynomial(fparsgt60[::-1])

sel = ((catalog['temperature_chi2'][ok] > 150) &
       (catalog['temperature_chi2'][ok] < 355))
fparsgt150 = np.polyfit(catalog['ratio303321'][ok][sel],
                        catalog['temperature_chi2'][ok][sel], 2)
p3 = np.polynomial.polynomial.Polynomial(fparsgt150[::-1])

root1 = np.polynomial.polynomial.polyroots((p1-p2).coef)
root1 = root1[(root1>0) & (root1<0.6)][0]
root2 = np.polynomial.polynomial.polyroots((p2-p3).coef)
root2 = root2[(root2>0) & (root2<0.6)][0]
#func = PiecewisePolynomial([0, root1, root2, 0.6],
#                           [p1.coef, p1.coef, p2.coef, p3.coef],)
func = lambda x: np.piecewise(x, [x<root1, (x>root1)&(x<root2), x>root2],
                              [p1, p2, p3])
log.info(" < {0}: {1}".format(root1, p1.coef))
log.info(" [{0}, {1}]: {2}".format(root1, root2, p2.coef))
log.info(" > {0}: {1}".format(root2, p3.coef))

fit_table = Table(
    [Column(name='MinBound', data=[0,root1,root2]),
     Column(name='MaxBound', data=[root1,root2,0.6]),
     Column(name='const',    data=[p.coef[0] if len(p.coef)>= 1 else 0
                                   for p in (p1,p2,p3)]),
     Column(name='xcoef',    data=[p.coef[1] if len(p.coef)>= 2 else 0
                                   for p in (p1,p2,p3)]),
     Column(name='x2coef',   data=[p.coef[2] if len(p.coef)>= 3 else 0
                                   for p in (p1,p2,p3)]),
    ])
fit_table.write(apath('piecewise_tvsratio_fit.ipac'), format='ascii.ipac')

x = np.linspace(0,0.6,100)
pl.plot(x, func(x), 'k--')
#ax2.plot(x, p1(x), 'k--')
#ax2.plot(x, p2(x), 'k--')
#ax2.plot(x, p3(x), 'k--')
ax2.set_xlim(0.,0.5)
ax2.set_ylim(10.,350)

fig3, ax3 = pl.subplots(num=3)
ax3.hist(sn[sn==sn], bins=50)
ax3.set_xlabel("Signal/Noise")

fig4, ax4 = pl.subplots(num=4)
for mask,color,alpha in masks_colors:
    ax4.errorbar(catalog['density_chi2'][mask], catalog['temperature_chi2'][mask],
                 #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                 #xerr=[catalog['elo_d'][mask], catalog['ehi_d'][mask]],
                 linestyle='none', capsize=0, alpha=alpha, marker='.',
                 color=color, linewidth=0.1)
    ax4.set_xlabel("Density")
    ax4.set_ylabel("Temperature")

fig5, ax5 = pl.subplots(num=5)
for mask,color,alpha in masks_colors:
    lon = catalog['x_cen'][mask]
    lon[lon>180] -= 360
    ax5.errorbar(lon, catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax5.set_xlabel("Galactic Longitude")
    ax5.set_ylabel("Temperature")

fig6, ax6 = pl.subplots(num=6)
for mask,color,alpha in masks_colors:
    ax6.errorbar(catalog['v_cen'][mask], catalog['temperature_chi2'][mask],
                #yerr=[catalog['elo_t'][mask], catalog['ehi_t'][mask]],
                linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
    ax6.set_xlabel("Centroid Velocity")
    ax6.set_ylabel("Temperature")

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


for ii in range(1,13):
    pl.figure(ii)
    if ii != 11:
        ax = pl.gca()
        ax.set_ylim(10, 125)
    pl.draw()

pl.draw()
pl.show()
