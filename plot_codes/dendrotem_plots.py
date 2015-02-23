import numpy as np
import pylab as pl
import matplotlib

from astropy.table import Table, Column
from astropy import log
from astropy import units as u
from scipy.interpolate import PiecewisePolynomial

from paths import hpath, apath, fpath, pcpath
from temperature_mapper import ph2cogrid, TemperatureMapper
from dendrograms import (catalog, catalog_sm, dend, dendsm)
import heating

matplotlib.rc_file(pcpath('pubfiguresrc'))

if 'tm' not in locals():
    tm = TemperatureMapper(logdensities=[3,4,5])
    tm2 = TemperatureMapper(logdensities=[4], deltav=20.0)
    tm3 = TemperatureMapper(logdensities=[4], deltav=1.0)

segmentdata = {'alpha': [(0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
               'blue': [(0.0, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
               'green': [(0.0, 0.0, 0.0), (0.5, 0.75, 0.75), (1.0, 0.0, 0.0)],
               'red': [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 1.0, 1.0)],
               'alpha':[(0.0,0.5,0.5), (0.5, 0.5, 0.5), (1.0, 0.5, 0.5)],
              }
cmap_rainbowvelo = matplotlib.colors.LinearSegmentedColormap(name='rgb',
                                                             segmentdata=segmentdata)

cbvmin,cbvmax = -50, 120
cblabel = r"$v_{LSR}$ (km s$^{-1}$)"


zipped = zip((catalog,catalog_sm,),#catalog321,catalog321_sm),
             (dend,dendsm,),#dend321,dend321sm),
             ('','_smooth',))#'_321ssel','_321sel_smooth'))


for cat,dendro,smooth in zipped[:1]:
    for ii in range(1,14):
        pl.figure(ii)
        pl.clf()

    sn = (cat['ratio303321']/cat['eratio303321'])
    sngt50 = sn > 50
    sn25_50 = (sn > 25) & (sn < 50)
    ok = (np.isfinite(sn) & (cat['Stot321'] < cat['Stot303']) & ~(cat['bad'] ==
                                                                  'True') &
          (cat['Smean321'] > 0) &
          (cat['e321'] > 0) &
          (~cat['IsNotH2CO']) & (~cat['IsAbsorption']))
    gt5 = (sn>5)

    hot = cat['temperature_chi2'] > 150

    is_leaf = np.array(cat['is_leaf'] == 'True')
    

    for ii in range(1,13): pl.figure(ii).clf()

    masks = (gt5 & ~sngt50 & ~sn25_50 & ok,
             sn25_50 & gt5 & ok,
             sngt50 & gt5 & ok,
             ok & ~gt5)
    leaf_masks = [mm for mask in masks for mm in (mask & is_leaf, mask & ~is_leaf)]
    # mask1 & leaf, mask1 & not leaf, mask2 & leaf, mask2 & not leaf....
    # Make the not-leaves be half as bright
    masks_colors = zip(leaf_masks,
                       ('b','b','g','g','r','r',    'k','k'),
                       (0.5,0.2, 0.6,0.3, 0.7,0.35, 0.3,0.15),
                       (8,7,9,8,10,9,5,4),
                      )

    fig1, ax1 = pl.subplots(num=1)
    for mask,color,alpha,markersize in masks_colors:
        ax1.errorbar(cat['area_exact'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax1.set_xscale('log')
        ax1.set_xlabel("Size Scale (pixels?)")
        ax1.set_ylabel("Temperature (K)")
    fig1.savefig(fpath('dendrotem/area_vs_temperature{0}.pdf'.format(smooth)))

    fig2, ax2 = pl.subplots(num=2)
    for mask,color,alpha,markersize in masks_colors:
        ax2.errorbar(cat['ratio303321'][mask], cat['temperature_chi2'][mask],
                     #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                     #xerr=[cat['eratio303321'][mask], cat['eratio303321'][mask]],
                     linestyle='none', capsize=0, alpha=alpha, marker='.',
                     markersize=10 if any(mask & is_leaf) else 5,
                     color=color, linewidth=0.3)
        ax2.set_xlabel("Ratio $S(3_{2,1}-2_{2,0})/S(3_{0,3}-2_{0,2})$")
        ax2.set_ylabel("Temperature (K)")
    fig2.savefig(fpath('dendrotem/ratio_vs_temperature{0}.pdf'.format(smooth)))

    if not hasattr(tm, 'temperatures'):
        tm.init()
    if not hasattr(tm2, 'temperatures'):
        tm2.init()
    if not hasattr(tm3, 'temperatures'):
        tm3.init()
    L1, = ax2.plot(tm.Xarr[1.2e-9]['ratio1'][1e3], tm.temperatures, 'k-.', label=r'$n(H_2)=10^3$ cm$^{-3}$', zorder=-5)
    L2, = ax2.plot(tm.Xarr[1.2e-9]['ratio1'][1e4], tm.temperatures, 'k--', label=r'$n(H_2)=10^4$ cm$^{-3}$', zorder=-5)
    L3, = ax2.plot(tm.Xarr[1.2e-9]['ratio1'][1e5], tm.temperatures, 'k:',  label=r'$n(H_2)=10^5$ cm$^{-3}$', zorder=-5)
    leg = pl.legend(loc='best')
    ax2.axis([0,0.55,10,200])
    fig2.savefig(fpath('dendrotem/ratio_vs_temperature{0}_modeloverlay.pdf'.format(smooth)))

    L4, = ax2.plot(tm2.Xarr[1.2e-9]['ratio1'][1e4],
                   tm2.temperatures, 'b--', alpha=0.2,
                   label=r'$n(H_2)=10^4$ cm$^{-3}$, $dv=20$ km s$^{-1}$', zorder=-10)
    L5, = ax2.plot(tm3.Xarr[1.2e-9]['ratio1'][1e4],
                   tm3.temperatures, 'r--', alpha=0.2,
                   label=r'$n(H_2)=10^4$ cm$^{-3}$, $dv=1$ km s$^{-1}$', zorder=-10)
    fig2.savefig(fpath('dendrotem/ratio_vs_temperature{0}_modeloverlay_dv.pdf'.format(smooth)))
    L4.set_visible(False)
    L5.set_visible(False)

    ax2.set_yscale('log')
    fig2.savefig(fpath('dendrotem/ratio_vs_temperature{0}_modeloverlay_log.pdf'.format(smooth)))
    ax2.set_yscale('linear')

    L1.set_visible(False)
    L2.set_visible(False)
    L3.set_visible(False)

        
    if cat is catalog:
        ## Determine approximate best-fit
        #sel = cat['temperature_chi2'][ok] < 70
        #fparslt60 = np.polyfit(cat['ratio303321'][ok][sel],
        #                       cat['temperature_chi2'][ok][sel], 2)
        #p1 = np.polynomial.polynomial.Polynomial(fparslt60[::-1])

        #sel = ((cat['temperature_chi2'][ok] > 70) &
        #       (cat['temperature_chi2'][ok] < 120) &
        #       (cat['ratio303321'][ok] < 0.6))
        #fparsgt60 = np.polyfit(cat['ratio303321'][ok][sel],
        #                       cat['temperature_chi2'][ok][sel], 1)
        #p2 = np.polynomial.polynomial.Polynomial(fparsgt60[::-1])

        #sel = ((cat['temperature_chi2'][ok] > 120) &
        #       (cat['temperature_chi2'][ok] < 355) &
        #       (cat['ratio303321'][ok] < 0.6))
        #fparsgt150 = np.polyfit(cat['ratio303321'][ok][sel],
        #                        cat['temperature_chi2'][ok][sel], 1)
        #p3 = np.polynomial.polynomial.Polynomial(fparsgt150[::-1])

        #root1 = np.polynomial.polynomial.polyroots((p1-p2).coef)
        #root1 = root1[(root1>0) & (root1<0.6)][0]
        #root2 = np.polynomial.polynomial.polyroots((p2-p3).coef)
        #root2 = root2[(root2>0) & (root2<0.6)][0]
        ##func = PiecewisePolynomial([0, root1, root2, 0.6],
        ##                           [p1.coef, p1.coef, p2.coef, p3.coef],)
        #func = lambda x: np.piecewise(x, [x<root1, (x>root1)&(x<root2), x>root2],
        #                              [p1, p2, p3])
        #log.info(" < {0}: {1}".format(root1, p1.coef))
        #log.info(" [{0}, {1}]: {2}".format(root1, root2, p2.coef))
        #log.info(" > {0}: {1}".format(root2, p3.coef))
        ## debug func = lambda x:x

        #fit_table = Table(
        #    [Column(name='MinBound', data=[0,root1,root2]),
        #     Column(name='MaxBound', data=[root1,root2,0.6]),
        #     Column(name='const',    data=[p.coef[0] if len(p.coef)>= 1 else 0
        #                                   for p in (p1,p2,p3)]),
        #     Column(name='xcoef',    data=[p.coef[1] if len(p.coef)>= 2 else 0
        #                                   for p in (p1,p2,p3)]),
        #     Column(name='x2coef',   data=[p.coef[2] if len(p.coef)>= 3 else 0
        #                                   for p in (p1,p2,p3)]),
        #    ])
        #fit_table.write(apath('piecewise_tvsratio_fit.ipac'), format='ascii.ipac')

        x = np.linspace(0,0.6,100)
        #l0, = ax2.plot(x, func(x), 'k-', alpha=0.5, zorder=-10)
        #ax2.set_xlim(0.,0.6)
        #ax2.set_ylim(0.,350)
        #fig2.savefig(fpath('dendrotem/ratio_vs_temperature_piecewise{0}.pdf'.format(smooth)))
        #l1, = ax2.plot(x, p1(x), 'k--', alpha=0.5)
        #l2, = ax2.plot(x, p2(x), 'k:', alpha=0.5)
        #l3, = ax2.plot(x, p3(x), 'k-.', alpha=0.5)
        #fig2.savefig(fpath('dendrotem/ratio_vs_temperature_piecewise_pieces{0}.pdf'.format(smooth)))
        #for ll in l1,l2,l3,l0:
        #    ll.set_visible(False)

        sel = cat['ratio303321'][ok] < 0.6
        pars = np.polyfit(cat['ratio303321'][ok][sel],
                          cat['temperature_chi2'][ok][sel],2)
        log.info("Polynomial parameters: {0}".format(pars))
        ax2.plot(x, np.polyval(pars, x), 'r--', alpha=0.5)
        ax2.axis([0,0.55,0,200])
        fig2.savefig(fpath('dendrotem/ratio_vs_temperature_powerlaw.pdf'))

    elif cat is catalog_sm:
        sel = cat['ratio303321'][ok] < 0.6
        pars = np.polyfit(cat['ratio303321'][ok][sel],
                          cat['temperature_chi2'][ok][sel],2)
        log.info("Polynomial parameters (smooth): {0}".format(pars))
        L, = ax2.plot(x, np.polyval(pars, x), 'k--', alpha=0.5,
                      label=r"$y={0:0.2f}x^2 {1:+0.2f}x {2:+0.2f}$".format(*pars))
        ax2.axis([0,0.55,0,200])
        pl.legend(loc='lower right')
        fig2.savefig(fpath('dendrotem/ratio_vs_temperature_powerlaw_smooth.pdf'))
        

    fig3, ax3 = pl.subplots(num=3)
    ax3.hist(sn[sn==sn], bins=50)
    ax3.set_xlabel("Signal/Noise")

    pl.figure(4).clf()
    fig4, (ax4a,ax4b) = pl.subplots(nrows=2,ncols=1,num=4)
    for mask,color,alpha,markersize in masks_colors:
        ax4a.errorbar(cat['density_chi2'][mask], cat['temperature_chi2'][mask],
                     #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                     #xerr=[cat['elo_d'][mask], cat['ehi_d'][mask]],
                     linestyle='none', capsize=0, alpha=alpha, marker='.',
                     color=color, linewidth=0.1)
        ax4a.set_xlabel("Density")
        ax4a.set_ylabel("Temperature")
        ax4a.set_xlim(3.0, 6.5)
        ax4b.errorbar(np.log10(cat['dustmindens'][mask]), cat['temperature_chi2'][mask],
                     #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                     #xerr=[cat['elo_d'][mask], cat['ehi_d'][mask]],
                     linestyle='none', capsize=0, alpha=alpha, marker='.',
                     color=color, linewidth=0.1)
        ax4b.set_xlabel("Density from dust")
        ax4b.set_ylabel("Temperature")
        ax4b.set_xlim(*ax4a.get_xlim())
    fig4.savefig(fpath("dendrotem/temperature_vs_density{0}.pdf".format(smooth)))

    fig5, ax5 = pl.subplots(num=5)
    lon = cat['x_cen']
    lon[lon>180] -= 360
    ax5.scatter(lon[hot], [149]*hot.sum(),
                c='r', alpha=0.3,
                s=1000*cat['Smean303'][hot],
                edgecolor='none',
                marker='^',)
    for mask,color,alpha,markersize in masks_colors:
        lon = cat['x_cen'][mask]
        lon[lon>180] -= 360
        ax5.scatter(lon, cat['temperature_chi2'][mask],
                    s=1000*cat['Smean303'][mask],
                    c=color,
                    edgecolor='none', alpha=alpha, marker='.')
        ax5.set_xlabel("Galactic Longitude")
        ax5.set_ylabel("Temperature (K)")
    ax5.set_ylim([0,150])
    fig5.savefig(fpath('dendrotem/temperature_vs_longitude{0}.pdf'.format(smooth)))

    fig5 = pl.figure(5)
    fig5.clf()
    ax5 = fig5.gca()
    lon = cat['x_cen']
    lon[lon>180] -= 360
    ax5.scatter(lon[hot], [149]*hot.sum(),
                c='r', alpha=0.3,
                s=1000*cat['Smean303'][hot],
                edgecolor='none',
                marker='^',)
    lon = cat['x_cen']
    lon[lon>180] -= 360
    vcen = cat['v_cen']/1e3
    color = cmap_rainbowvelo((vcen-cbvmin)/(cbvmax-cbvmin))
    sc = ax5.scatter(lon[is_leaf&ok], cat['temperature_chi2'][is_leaf&ok],
                     s=1000*cat['Smean303'][is_leaf&ok],
                     c=color[is_leaf&ok],
                     edgecolor='none', marker='.')
    sc2 = ax5.scatter(lon[(~is_leaf)&ok], cat['temperature_chi2'][(~is_leaf)&ok],
                      s=1000*cat['Smean303'][(~is_leaf)&ok],
                      c=color[(~is_leaf)&ok],
                      edgecolor='none', alpha=0.2, marker='.')
    ax5.set_xlabel("Galactic Longitude")
    ax5.set_ylabel("Temperature (K)")
    ax5.set_ylim([0,150])
    ax5.set_xlim([1.7,-0.6])
    sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=cbvmin, vmax=cbvmax),
                                      cmap=cmap_rainbowvelo)
    sm._A = []
    cb = fig5.colorbar(sm)
    cb.set_label(cblabel)
    fig5.savefig(fpath('dendrotem/temperature_vs_longitude_velocolor{0}.pdf'.format(smooth)))


    fig6, ax6 = pl.subplots(num=6)
    for mask,color,alpha,markersize in masks_colors:
        ax6.errorbar(cat['v_cen'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax6.set_xlabel("Centroid Velocity")
        ax6.set_ylabel("Temperature (K)")

    fig7, ax7 = pl.subplots(num=7)
    for mask,color,alpha,markersize in masks_colors:
        ax7.errorbar(cat['radius'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax7.set_xscale('log')
        ax7.set_xlabel("Radius (pixels?)")
        ax7.set_ylabel("Temperature (K)")

    fig8, ax8 = pl.subplots(num=8)
    for mask,color,alpha,markersize in masks_colors:
        ax8.errorbar(cat['Smean303'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    #xerr=[cat['e303'][mask], cat['e303'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax8.set_xlabel("$S(3_{0,3}-2_{0,2})$ (K)")
        ax8.set_ylabel("Temperature (K)")

    fig9, ax9 = pl.subplots(num=9)
    for mask,color,alpha,markersize in masks_colors:
        ax9.errorbar(cat['Smean321'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    #xerr=[cat['e321'][mask], cat['e321'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax9.set_xlabel("$S(3_{2,1}-2_{2,0})$ (K)")
        ax9.set_ylabel("Temperature (K)")

    fig10, ax10 = pl.subplots(num=10)
    for mask,color,alpha,markersize in masks_colors:
        ax10.errorbar(cat['13comean'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax10.set_xlabel("${S}(^{13}$CO) (K)")
        ax10.set_ylabel("Temperature (K)")

    fig11, ax11 = pl.subplots(num=11)
    for mask,color,alpha,markersize in masks_colors:
        ax11.errorbar(cat['13comean'][mask], cat['Smean303'][mask],
                    #yerr=[cat['e303'][mask], cat['e303'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax11.set_xlabel("$\\bar{S}(^{13}$CO) (K)")
        ax11.set_ylabel("$S(3_{0,3}-2_{0,2})$ (K)")

    fig12, ax12 = pl.subplots(num=12)
    ax12.errorbar(cat['v_rms'][hot]*np.sqrt(8*np.log(2)), [149]*hot.sum(),
                  lolims=True, linestyle='none', capsize=0, alpha=0.3,
                  marker='^', color='r')
    for mask,color,alpha,markersize in masks_colors:
        ax12.errorbar(cat['v_rms'][mask]*np.sqrt(8*np.log(2)), cat['temperature_chi2'][mask],
                      #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                      markersize=10 if any(mask & is_leaf) else 5,
                      markeredgecolor='none',
                      linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax12.set_xlabel(r"Line FWHM (km s$^{-1}$)")
        ax12.set_ylabel("Temperature (K)")

    linewidths = np.linspace(0,cat['v_rms'].max())*u.km/u.s
    ax12.plot(linewidths*2.35, [heating.tkin_all(10**4*u.cm**-3, sigma, 10*u.pc,
                                                5*u.km/u.s/u.pc, 30*u.K)
                               for sigma in linewidths],
            linestyle='--', color='k', label='$n=10^4$ cm$^{-3}$', zorder=-5)
    ax12.plot(linewidths*2.35, [heating.tkin_all(10**4*u.cm**-3, sigma, 10*u.pc,
                                                1*u.km/u.s/u.pc, 30*u.K)
                               for sigma in linewidths],
            linestyle='--', color='r', label='$n=10^4$ cm$^{-3}$, $dv/dr=1$', zorder=-5, linewidth=2, alpha=0.5)
    ax12.plot(linewidths*2.35, [heating.tkin_all(10**4*u.cm**-3, sigma, 20*u.pc,
                                                5*u.km/u.s/u.pc, 30*u.K)
                               for sigma in linewidths],
             linestyle='--', color='b', label='$n=10^4$ cm$^{-3}$, $L=20$ pc', zorder=-5, alpha=0.5, linewidth=2)
    ax12.plot(linewidths*2.35, [heating.tkin_all(10**5*u.cm**-3, sigma, 10*u.pc,
                                                5*u.km/u.s/u.pc, 30*u.K)
                               for sigma in linewidths],
             linestyle=':', color='k', label='$n=10^5$ cm$^{-3}$', zorder=-5)
    ax12.plot(linewidths*2.35, [heating.tkin_all(10**6*u.cm**-3, sigma, 10*u.pc,
                                                5*u.km/u.s/u.pc, 30*u.K)
                               for sigma in linewidths],
            linestyle='-.', color='k', label='$n=10^6$ cm$^{-3}$', zorder=-5)
    ax12.plot(linewidths*2.35, [heating.tkin_all(10**5*u.cm**-3, sigma, 10*u.pc,
                                                5*u.km/u.s/u.pc, 30*u.K, crir=1e-15*u.s**-1)
                               for sigma in linewidths],
             linestyle='-', color='g', label='$n=10^5$ cm$^{-3}$, $\zeta_{CR}=10^{-15}$ s$^{-1}$', zorder=-10, alpha=0.25, linewidth=4)
    ax12.plot(linewidths*2.35, [heating.tkin_all(10**5*u.cm**-3, sigma, 10*u.pc,
                                                5*u.km/u.s/u.pc, 30*u.K, crir=1e-14*u.s**-1)
                               for sigma in linewidths],
              linestyle=':', color='purple', label='$n=10^5$ cm$^{-3}$, $\zeta_{CR}=10^{-14}$ s$^{-1}$', zorder=-10, alpha=0.25, linewidth=4)

    ax12.set_ylim([0,150])
    fig12.savefig(fpath('dendrotem/temperature_vs_rmsvelocity{0}.pdf'.format(smooth)))
    wide = cat['v_rms'] > 20/np.sqrt(8*np.log(2))
    ax12.errorbar([24.5] * (wide & is_leaf).sum(),
                  cat['temperature_chi2'][wide&is_leaf],
                  lolims=True, linestyle='none', capsize=0, alpha=0.3,
                  markersize=10,
                  marker='>', color='r')
    ax12.errorbar([24.5] * (wide & ~is_leaf).sum(),
                  cat['temperature_chi2'][wide&(~is_leaf)],
                  lolims=True, linestyle='none', capsize=0, alpha=0.1,
                  markersize=5,
                  marker='>', color='r')
    ax12.set_xlim([0,25])
    fig12.savefig(fpath('dendrotem/temperature_vs_rmsvelocity_xzoom{0}.pdf'.format(smooth)))

    fig22 = pl.figure(22)
    fig22.clf()
    ax22 = fig22.gca()
    ax22.errorbar(cat['higaldusttem'][hot], cat['tmin1sig_chi2'][hot],
                  lolims=True, linestyle='none', capsize=0, alpha=alpha,
                  markeredgecolor='none',
                  marker='^', color='r')
    for mask,color,alpha,markersize in masks_colors:
        ax22.errorbar(cat['higaldusttem'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                      markersize=markersize*(2 if any(mask & is_leaf) else 1),
                      linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        mask = mask & (cat['tmin1sig_chi2'] < cat['higaldusttem'])
        ax22.plot(cat['higaldusttem'][mask], cat['temperature_chi2'][mask],
                  markersize=markersize*(2 if any(mask & is_leaf) else 1)*1.1,
                  linestyle='none',
                  alpha=alpha,
                  marker='o',
                  markerfacecolor='none',
                  markeredgecolor=color)
    ax22.plot([14,38], [14,38], 'k--', zorder=-5, alpha=0.5)
    ax22.set_xlim([13,35])
    ax22.set_ylim([13,150])
    ax22.set_xlabel("HiGal Dust Temperature (K)")
    ax22.set_ylabel("Temperature (K)")
    fig22.savefig(fpath('dendrotem/temperature_vs_dusttem{0}.pdf'.format(smooth)))

    fig24 = pl.figure(24)
    fig24.clf()
    ax24 = fig24.gca()
    hot_lolim = cat['tmin1sig_chi2'] > 150
    ax24.errorbar(10**cat['logh2column'][hot|hot_lolim], [149]*(hot|hot_lolim).sum(),
                  lolims=True, linestyle='none', capsize=0, alpha=alpha,
                  marker='^', color='r')
    for mask,color,alpha,markersize in masks_colors:
        ax24.errorbar(10**cat['logh2column'][mask&~hot_lolim],
                      cat['temperature_chi2'][mask&~hot_lolim],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax24.plot([14,38], [14,38], 'k--')
        ax24.set_ylim([13,150])
        ax24.set_xlabel("HiGal Dust Column Density")
        ax24.set_ylabel("Temperature (K)")
    fig24.savefig(fpath('dendrotem/temperature_vs_dustcol{0}.pdf'.format(smooth)))

    fig26 = pl.figure(26)
    fig26.clf()
    ax26 = fig26.gca()
    mask = (hot|hot_lolim) & is_leaf
    ax26.errorbar(cat['tkin_turb'][mask],
                  cat['elo_t'][mask],
                  lolims=True, linestyle='none', capsize=0, alpha=0.3,
                  marker='^', color='r')
    for mask,color,alpha,markersize in masks_colors:
        mask = (mask & (~hot_lolim) & is_leaf)
        if mask.sum() == 0: continue
        ax26.errorbar(cat['tkin_turb'][mask],
                      cat['temperature_chi2'][mask],
                      yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                      linestyle='none', capsize=0, alpha=0.5*alpha, marker='',
                      linewidth=0.5,
                      color=color, markeredgecolor='none')
        ax26.plot(cat['tkin_turb'][mask],
                  cat['temperature_chi2'][mask],
                  linestyle='none',  alpha=alpha, marker='o', color=color,
                  markeredgecolor='none')

        # Highlight those inconsistent with the curve
        mask = mask & (cat['tmin1sig_chi2'] > cat['tkin_turb'])
        ax26.plot(cat['tkin_turb'][mask],
                  cat['temperature_chi2'][mask],
                  linestyle='none',  alpha=alpha, marker='o',
                  markersize=10,
                  markerfacecolor='none',
                  markeredgewidth=0.3,
                  markeredgecolor=color)

    ax26.plot([0,200], [0,200], 'k--', alpha=0.5, zorder=-5)
    ax26.set_ylim([0,155])
    ax26.set_xlim([cat['tkin_turb'][is_leaf].min()-2,cat['tkin_turb'][is_leaf].max()+2])
    ax26.set_xlabel("Turbulence-driven Temperature (K)")
    ax26.set_ylabel("H$_2$CO Temperature (K)")

    fig26.savefig(fpath('dendrotem/temperature_vs_turbtemperature{0}.pdf'.format(smooth)))

    for mask,color,alpha,markersize in masks_colors:
        mask = (mask & (~hot_lolim) & ~is_leaf)
        if mask.sum() == 0: continue
        ax26.errorbar(cat['tkin_turb'][mask],
                      cat['temperature_chi2'][mask],
                      yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                      linestyle='none', capsize=0, alpha=0.5*alpha, marker='',
                      linewidth=0.5,
                      color=color, markeredgecolor='none')
        ax26.plot(cat['tkin_turb'][mask],
                  cat['temperature_chi2'][mask],
                  linestyle='none',  alpha=alpha, marker='o', color=color,
                  markeredgecolor='none')

        # Highlight those inconsistent with the curve
        mask = mask & (cat['tmin1sig_chi2'] > cat['tkin_turb'])
        ax26.plot(cat['tkin_turb'][mask],
                  cat['temperature_chi2'][mask],
                  linestyle='none',  alpha=alpha, marker='o',
                  markersize=10,
                  markerfacecolor='none',
                  markeredgewidth=0.3,
                  markeredgecolor=color)

    ax26.set_ylim([0,155])
    ax26.set_xlim([cat['tkin_turb'][is_leaf].min()-2,cat['tkin_turb'][is_leaf].max()+2])

    fig26.savefig(fpath('dendrotem/temperature_vs_turbtemperature_withtrunks{0}.pdf'.format(smooth)))

    fig25 = pl.figure(25)
    fig25.clf()
    ax25 = fig25.gca()
    tem_to_color = pl.cm.RdYlBu_r((cat['temperature_chi2']-15)/(200-15.))
    sc = ax25.scatter(cat['ratio303321'][is_leaf], np.log10(cat['dustmindens'][is_leaf]),
                      c=cat['temperature_chi2'][is_leaf], cmap=pl.cm.RdYlBu_r,
                      marker='o', vmin=15, vmax=200, alpha=0.8,
                      lw=0, s=40, edgecolors='k')
    sc2 = ax25.scatter(cat['ratio303321'][~is_leaf], np.log10(cat['dustmindens'][~is_leaf]),
                       c=cat['temperature_chi2'][~is_leaf], cmap=pl.cm.RdYlBu_r,
                       marker='o', vmin=15, vmax=200, alpha=0.2,
                       lw=0, s=25, edgecolors='k')
    ax25.set_xlim(0, 0.6)
    ax25.set_ylabel("HiGal Dust-derived Density")
    ax25.set_xlabel("Ratio $R_1$")
    #ax25.set_axis_bgcolor((0.6,0.6,0.6))
    pl.colorbar(sc)
    fig25.savefig(fpath('dendrotem/ratio_vs_density{0}.pdf'.format(smooth)))


    fig13, ax13 = pl.subplots(num=13)
    lon=cat['x_cen']
    lon[lon>180] -= 360
    pts13 = ax13.scatter(lon, cat['y_cen'],
                         c=cat['temperature_chi2'],
                         alpha=0.3, marker='o',
                         norm=matplotlib.colors.PowerNorm(0.5),
                         cmap=matplotlib.cm.RdYlBu_r,
                         s=(cat['area_exact'])**0.5)
    ax13.clear()
    pn = matplotlib.colors.PowerNorm(0.5, vmax=np.nanmax(cat['temperature_chi2']), vmin=np.nanmin(cat['temperature_chi2']))
    y = pn(cat['temperature_chi2'])
    colors = matplotlib.cm.RdYlBu_r(y)
    colors[:,3] = np.nan_to_num(np.array((pn(cat['temperature_chi2'])+1)/2.))
    sz = (cat['area_exact'])**0.5
    colors[:,3][sz>100] /= 3
    ln = matplotlib.colors.LogNorm(vmax=100, vmin=3, clip=True)
    colors[:,3] = np.nan_to_num(ln(sn))
    pts13b = ax13.scatter(lon, cat['y_cen'], c=colors, marker='o',
                          edgecolor='none',
                          norm=matplotlib.colors.PowerNorm(0.5), alpha=0.75,
                          s=sz)
    ax13.set_xlabel("Galactic Longitude")
    ax13.set_ylabel("Galactic Latitude")
    cb13 = pl.colorbar(pts13)
    cb13.set_label(r'Temperature')
    ax13.set_xlim(1.7,-0.6)
    ax13.set_axis_bgcolor((0.6,0.6,0.6))



    brick_coords = [(262-64)/(2 if 'smooth' in smooth else 1),143,725]
    sgra_coords  = [(239-64)/(2 if 'smooth' in smooth else 1),97,907]
    sgrb2_coords = [(278-64)/(2 if 'smooth' in smooth else 1),117,522]

    try:
        brick = dendro.structure_at(brick_coords).ancestor
        sgra = dendro.structure_at(sgra_coords).ancestor
        sgrb2 = dendro.structure_at(sgrb2_coords).ancestor
    except AttributeError:
        log.warning("Failed to find one of the sources.")
        continue
    brick_leaves = [obj for obj in brick.descendants if obj.is_leaf]
    sgra_leaves = [obj for obj in sgra.descendants if obj.is_leaf]
    sgrb2_leaves = [obj for obj in sgrb2.descendants if obj.is_leaf]

    fig14 = pl.figure(14)
    fig14.clf()
    ax14 = fig14.gca()
    def dendroplot(axis=ax14, axname1='area_exact', axname2='ratio303321',
                   axscale1=1.,
                   axscale2=1.,
                   leaves_list=[sgra_leaves],
                   # r, b, g
                   color_list=['#CC4444', '#4444CC', '#44CC44'],
                   highlight_monotonic=True,
                   marker='s',
                   marker2=None,
                   linestyle='-', **kwargs):
        for leaves, color in zip(leaves_list,color_list):
            for leaf in leaves:
                xax,yax = ([cat[leaf.idx][axname1]*axscale1],
                           [cat[leaf.idx][axname2]*axscale2])
                axis.plot(xax, yax, marker, color=color, markeredgecolor='none', alpha=0.5)
                obj = leaf.parent
                while obj.parent:
                    xax.append(cat[obj.idx][axname1]*axscale1)
                    yax.append(cat[obj.idx][axname2]*axscale2)
                    obj = obj.parent
                if np.any(np.isnan(yax)):
                    ok = ~np.isnan(yax)
                    axis.plot(np.array(xax)[ok], np.array(yax)[ok], alpha=0.5,
                              label=leaf.idx, color='b', zorder=5,
                              linestyle=linestyle, marker=marker2, **kwargs)
                else:
                    axis.plot(xax, yax, alpha=0.1, label=leaf.idx, color=color,
                              zorder=5, linestyle=linestyle, marker=marker2,
                              **kwargs)
                if highlight_monotonic:
                    signs = np.sign(np.diff(yax))
                    if np.all(signs==1) or np.all(signs==-1):
                        axis.plot(xax, yax, alpha=0.1, linewidth=5, zorder=0, color='g')
    dendroplot()
    ax14.set_xscale('log')
    ax14.set_xlabel("Area (square arcseconds)")
    ax14.set_ylabel("Ratio 321/303")
    fig14.savefig(fpath('dendrotem/sgra_ratio_vs_sizescale{0}.pdf'.format(smooth)))


    fig15 = pl.figure(15)
    fig15.clf()
    ax15 = fig15.gca()
    dendroplot(axis=ax15, axname2='Smean303')
    ax15.set_xscale('log')
    ax15.set_xlabel("Area (square arcseconds)")
    ax15.set_ylabel(r"$\bar{S_\nu}(3_{03}-2_{02})$")
    fig15.savefig(fpath('dendrotem/flux303mean_vs_area{0}.pdf'.format(smooth)))

    fig16 = pl.figure(16)
    fig16.clf()
    ax16 = fig16.gca()
    dendroplot(axis=ax16, axname2='Stot303')
    ax16.set_xscale('log')
    ax16.set_yscale('log')
    ax16.set_xlabel("Area (square arcseconds)")
    ax16.set_ylabel(r"$\Sigma S_\nu(3_{03}-2_{02})$")
    fig16.savefig(fpath('dendrotem/flux303sum_vs_area{0}.pdf'.format(smooth)))


    fig17 = pl.figure(17)
    fig17.clf()
    ax17 = fig17.gca()
    dendroplot(axis=ax17, axname2='temperature_chi2', leaves_list=[sgra_leaves,sgrb2_leaves])
    ax17.set_xscale('log')
    ax17.set_xlabel("Area (square arcseconds)")
    ax17.set_ylabel(r"Temperature (K)")
    fig17.savefig(fpath('dendrotem/temperature_vs_area{0}.pdf'.format(smooth)))

    fig17.clf()
    ax17 = fig17.gca()
    dendroplot(axis=ax17, axname1='Smean303', axname2='temperature_chi2', leaves_list=[sgra_leaves,sgrb2_leaves])
    ax17.set_xscale('log')
    ax17.set_xlabel(r"$\bar{S_\nu}(3_{03}-2_{02})$")
    ax17.set_ylabel(r"Temperature (K)")
    fig17.savefig(fpath('dendrotem/temperature_vs_flux{0}.pdf'.format(smooth)))

    fig18 = pl.figure(18)
    fig18.clf()
    ax18 = fig18.gca()
    dendroplot(leaves_list=[sgra_leaves, brick_leaves], axname2='temperature_chi2', axis=ax18)
    ax18.set_xscale('log')
    ax18.set_xlabel("Area (square arcseconds)")
    ax18.set_ylabel("Temperature (K)")
    fig18.savefig(fpath('dendrotem/all_temperature_vs_sizescale{0}.pdf'.format(smooth)))

    fig19 = pl.figure(19)
    fig19.clf()
    ax19 = fig19.gca()
    ax19.errorbar(cat['Smean303'], cat['Smean321'],
                    linestyle='none', capsize=0, alpha=0.2, marker='.', color='k')
    dendroplot(leaves_list=[sgra_leaves, brick_leaves, sgrb2_leaves], axname1='Smean303',
               axname2='Smean321', axis=ax19, highlight_monotonic=False, linestyle='none',
               marker='.', marker2='.')
    ax19.set_xlabel(r"$\bar{S}_\nu(3_{03}-2_{02})$")
    ax19.set_ylabel(r"$\bar{S}_\nu(3_{21}-2_{20})$")
    ax19.axis([0.1,2,0.01,1])
    ax19.set_xscale('log')
    ax19.set_yscale('log')
    fig19.savefig(fpath('dendrotem/S303vsS321{0}.pdf'.format(smooth)))

    fig20 = pl.figure(20)
    fig20.clf()
    ax20 = fig20.gca()
    dendroplot(leaves_list=[sgra_leaves, brick_leaves, sgrb2_leaves], axname1='area_exact',
               axname2='v_cen', axis=ax20, highlight_monotonic=False,)
    ax20.set_ylabel(r"$v_{cen}$")
    ax20.set_xlabel(r"Area")
    ax20.set_xscale('log')
    fig20.savefig(fpath('dendrotem/vlsr_vs_area{0}.pdf'.format(smooth)))


    fig21 = pl.figure(21)
    fig21.clf()
    ax21 = fig21.gca()
    dendroplot(leaves_list=[sgra_leaves, brick_leaves, sgrb2_leaves], axname1='area_exact',
               axname2='v_rms', axis=ax21, highlight_monotonic=False,)
    ax21.set_ylabel(r"$v_{rms}=\sigma_v$")
    ax21.set_xlabel(r"Area")
    ax21.set_xscale('log')
    fig21.savefig(fpath('dendrotem/vrms_vs_area{0}.pdf'.format(smooth)))

    fig23 = pl.figure(23)
    fig23.clf()
    ax23 = fig23.gca()

    ax23.errorbar(cat['v_rms'][hot]*np.sqrt(8*np.log(2)), [149]*hot.sum(),
                  lolims=True, linestyle='none', capsize=0, alpha=alpha,
                  marker='^', color='r')
    for mask,color,alpha,markersize in masks_colors:
        ax23.errorbar(cat['v_rms'][mask]*np.sqrt(8*np.log(2)), cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax23.set_xlabel(r"Line FWHM (km s$^{-1}$)")
        ax23.set_ylabel("H$_2$CO Temperature (K)")
    ax23.set_ylim([0,150])

    # r,b,g = sgra, brick, sgrb2
    dendroplot(leaves_list=[sgra_leaves, brick_leaves, sgrb2_leaves],
               axname1='v_rms', axname2='temperature_chi2', axis=ax23,
               axscale1=np.sqrt(8*np.log(2)),
               linestyle='none', marker='s', markerfacecolor='none',
               highlight_monotonic=False,)
    ax23.set_xlabel(r"Line FWHM (km s$^{-1}$)")
    ax23.set_ylabel(r"H$_2$CO Temperature (K)")
    fig23.savefig(fpath('dendrotem/dendro_temperature_vs_rmsvelocity{0}.pdf'.format(smooth)))




    for ii in range(1,13):
        pl.figure(ii)
        if ii not in (11,12,13,14,15):
            ax = pl.gca()
            ax.set_ylim(10, 125)
        pl.draw()

    pl.close(27)


    dview = dendro.viewer()
    structure = dendro.structure_at(brick_coords).ancestor
    dview.hub.select(1, structure)
    dview.ax_image.axis((710,740,131,154))
    dview.slice_slider.set_val(brick_coords[0])
    dview.fig.savefig(fpath('dendrotem/dendrogram_viewer_brick{0}.pdf'.format(smooth)))

    dview = dendro.viewer()
    structure = dendro.structure_at(sgra_coords).ancestor
    dview.hub.select(1, structure)
    dview.ax_image.axis([sgra_coords[2]-20, sgra_coords[2]+20, sgra_coords[1]-20, sgra_coords[1]+20])
    dview.slice_slider.set_val(sgra_coords[0])
    dview.fig.savefig(fpath('dendrotem/dendrogram_viewer_sgra{0}.pdf'.format(smooth)))


    # SgrA total, 20kms, 50kms
    #catalog_sm[np.array([126,247,346])]['_idx','ratio303321','dustmindens','temperature_chi2'].pprint()

    

    pl.draw()
    pl.show()


    # catalog_sm[np.abs(catalog_sm['temperature_chi2']-catalog_sm['higaldusttem'])/catalog_sm['higaldusttem'] < 1.5].pprint()

    def zoom_to(idx, dendro=dendro, viewer=dview, button=1):
        viewer.hub.select(button, dendro[idx])
        w = dendro.wcs
        sw = w.sub([3])
        cw = w.sub([1,2])
        row = cat[idx]
        viewer.slice_slider.set_val(sw.wcs_world2pix((row['v_cen'],),0)[0])
        edges = cw.wcs_world2pix([(row['x_cen'] - row['major_sigma']/3600.*3,
                                   row['y_cen'] - row['major_sigma']/3600.*3),
                                  (row['x_cen'] + row['major_sigma']/3600.*3,
                                   row['y_cen'] + row['major_sigma']/3600.*3)],
                                 0)
        x0,x1 = edges[:,0] if edges[1,0] > edges[0,0] else edges[::-1,0]
        y0,y1 = edges[:,1] if edges[1,1] > edges[0,1] else edges[::-1,1]
        viewer.ax_image.axis((x0,x1,y0,y1))
        viewer.ax_image.axis((x0,x1,y0,y1))
        pl.draw()
        pl.show()
