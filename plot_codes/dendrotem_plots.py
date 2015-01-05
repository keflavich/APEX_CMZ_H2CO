import numpy as np
import pylab as pl
import matplotlib

from paths import hpath, apath, fpath, pcpath
from astropy.table import Table, Column
from astropy import log
from scipy.interpolate import PiecewisePolynomial

from dendrograms import (catalog, catalog_sm, dend, dendsm, dend321, dend321sm,
                         catalog321, catalog321_sm)

matplotlib.rc_file(pcpath('pubfiguresrc'))

zipped = zip((catalog,catalog_sm,),#catalog321,catalog321_sm),
             (dend,dendsm,),#dend321,dend321sm),
             ('','_smooth',))#'_321ssel','_321sel_smooth'))

for cat,dendro,smooth in zipped:
    for ii in range(1,14):
        pl.figure(ii)
        pl.clf()

    sn = (cat['ratio303321']/cat['eratio303321'])
    sngt50 = sn > 50
    sn25_50 = (sn > 25) & (sn < 50)
    ok = np.isfinite(sn) & (cat['Stot321'] < cat['Stot303'])
    gt5 = (sn>5)
    

    for ii in range(1,13): pl.figure(ii).clf()

    masks_colors = zip((gt5 & ~sngt50 & ~sn25_50, sn25_50 & gt5, sngt50 & gt5, ok & ~gt5),
                       ('b','g','r','k'),
                       (0.2,0.3,0.4,0.1),
                      )

    fig1, ax1 = pl.subplots(num=1)
    for mask,color,alpha in masks_colors:
        ax1.errorbar(cat['area_exact'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax1.set_xscale('log')
        ax1.set_xlabel("Size Scale (pixels?)")
        ax1.set_ylabel("Temperature")
    fig1.savefig(fpath('dendrotem/area_vs_temperature{0}.pdf'.format(smooth)))

    fig2, ax2 = pl.subplots(num=2)
    for mask,color,alpha in masks_colors:
        ax2.errorbar(cat['ratio303321'][mask], cat['temperature_chi2'][mask],
                     #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                     #xerr=[cat['eratio303321'][mask], cat['eratio303321'][mask]],
                     linestyle='none', capsize=0, alpha=alpha, marker='.',
                     markersize=10,
                     color=color, linewidth=0.3)
        ax2.set_xlabel("Ratio $S(3_{2,1}-2_{2,0})/S(3_{0,3}-2_{0,2})$")
        ax2.set_ylabel("Temperature [K]")
    fig2.savefig(fpath('dendrotem/ratio_vs_temperature{0}.pdf'.format(smooth)))
        
    if cat is catalog:
        # Determine approximate best-fit
        sel = cat['temperature_chi2'][ok] < 70
        fparslt60 = np.polyfit(cat['ratio303321'][ok][sel],
                               cat['temperature_chi2'][ok][sel], 2)
        p1 = np.polynomial.polynomial.Polynomial(fparslt60[::-1])

        sel = ((cat['temperature_chi2'][ok] > 70) &
               (cat['temperature_chi2'][ok] < 120) &
               (cat['ratio303321'][ok] < 0.6))
        fparsgt60 = np.polyfit(cat['ratio303321'][ok][sel],
                               cat['temperature_chi2'][ok][sel], 1)
        p2 = np.polynomial.polynomial.Polynomial(fparsgt60[::-1])

        sel = ((cat['temperature_chi2'][ok] > 120) &
               (cat['temperature_chi2'][ok] < 355) &
               (cat['ratio303321'][ok] < 0.6))
        fparsgt150 = np.polyfit(cat['ratio303321'][ok][sel],
                                cat['temperature_chi2'][ok][sel], 1)
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
        ax2.plot(x, func(x), 'k-', alpha=0.5, zorder=-10)
        ax2.set_xlim(0.,0.6)
        ax2.set_ylim(0.,350)
        fig2.savefig(fpath('dendrotem/ratio_vs_temperature_piecewise{0}.pdf'.format(smooth)))
        ax2.plot(x, p1(x), 'k--', alpha=0.5)
        ax2.plot(x, p2(x), 'k:', alpha=0.5)
        ax2.plot(x, p3(x), 'k-.', alpha=0.5)
        fig2.savefig(fpath('dendrotem/ratio_vs_temperature_piecewise_pieces{0}.pdf'.format(smooth)))

    fig3, ax3 = pl.subplots(num=3)
    ax3.hist(sn[sn==sn], bins=50)
    ax3.set_xlabel("Signal/Noise")

    fig4, ax4 = pl.subplots(num=4)
    for mask,color,alpha in masks_colors:
        ax4.errorbar(cat['density_chi2'][mask], cat['temperature_chi2'][mask],
                     #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                     #xerr=[cat['elo_d'][mask], cat['ehi_d'][mask]],
                     linestyle='none', capsize=0, alpha=alpha, marker='.',
                     color=color, linewidth=0.1)
        ax4.set_xlabel("Density")
        ax4.set_ylabel("Temperature")

    fig5, ax5 = pl.subplots(num=5)
    for mask,color,alpha in masks_colors:
        lon = cat['x_cen'][mask]
        lon[lon>180] -= 360
        ax5.errorbar(lon, cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax5.set_xlabel("Galactic Longitude")
        ax5.set_ylabel("Temperature")
    fig5.savefig(fpath('dendrotem/temperature_vs_longitude{0}.pdf'.format(smooth)))

    fig6, ax6 = pl.subplots(num=6)
    for mask,color,alpha in masks_colors:
        ax6.errorbar(cat['v_cen'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax6.set_xlabel("Centroid Velocity")
        ax6.set_ylabel("Temperature")

    fig7, ax7 = pl.subplots(num=7)
    for mask,color,alpha in masks_colors:
        ax7.errorbar(cat['radius'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax7.set_xscale('log')
        ax7.set_xlabel("Radius (pixels?)")
        ax7.set_ylabel("Temperature")

    fig8, ax8 = pl.subplots(num=8)
    for mask,color,alpha in masks_colors:
        ax8.errorbar(cat['Smean303'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    #xerr=[cat['e303'][mask], cat['e303'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax8.set_xlabel("$S(3_{0,3}-2_{0,2})$ (K)")
        ax8.set_ylabel("Temperature [K]")

    fig9, ax9 = pl.subplots(num=9)
    for mask,color,alpha in masks_colors:
        ax9.errorbar(cat['Smean321'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    #xerr=[cat['e321'][mask], cat['e321'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax9.set_xlabel("$S(3_{2,1}-2_{2,0})$ (K)")
        ax9.set_ylabel("Temperature [K]")

    fig10, ax10 = pl.subplots(num=10)
    for mask,color,alpha in masks_colors:
        ax10.errorbar(cat['13comean'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax10.set_xlabel("${S}(^{13}$CO) (K)")
        ax10.set_ylabel("Temperature [K]")

    fig11, ax11 = pl.subplots(num=11)
    for mask,color,alpha in masks_colors:
        ax11.errorbar(cat['13comean'][mask], cat['Smean303'][mask],
                    #yerr=[cat['e303'][mask], cat['e303'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax11.set_xlabel("$\\bar{S}(^{13}$CO) (K)")
        ax11.set_ylabel("$S(3_{0,3}-2_{0,2})$ (K)")

    fig12, ax12 = pl.subplots(num=12)
    hot = cat['temperature_chi2'] > 150
    ax12.errorbar(cat['v_rms'][hot], [149]*hot.sum(),
                  lolims=True, linestyle='none', capsize=0, alpha=alpha,
                  marker='^', color='r')
    for mask,color,alpha in masks_colors:
        ax12.errorbar(cat['v_rms'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax12.set_xlabel("RMS Velocity")
        ax12.set_ylabel("Temperature [K]")
    ax12.set_ylim([0,150])
    fig12.savefig(fpath('dendrotem/temperature_vs_rmsvelocity{0}.pdf'.format(smooth)))

    fig22 = pl.figure(22)
    fig22.clf()
    ax22 = fig22.gca()
    hot = cat['temperature_chi2'] > 150
    ax22.errorbar(cat['higaldusttem'][hot], [149]*hot.sum(),
                  lolims=True, linestyle='none', capsize=0, alpha=alpha,
                  marker='^', color='r')
    for mask,color,alpha in masks_colors:
        ax22.errorbar(cat['higaldusttem'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax22.plot([14,38], [14,38], 'k--')
        ax22.set_xlim([13,38])
        ax22.set_ylim([13,150])
        ax22.set_xlabel("HiGal Dust Temperature [K]")
        ax22.set_ylabel("Temperature [K]")
    fig22.savefig(fpath('dendrotem/temperature_vs_dusttem{0}.pdf'.format(smooth)))

    fig13, ax13 = pl.subplots(num=13)
    lon=cat['x_cen']
    lon[lon>180] -= 360
    pts13 = ax13.scatter(lon, cat['y_cen'],
                         c=cat['temperature_chi2'],
                         alpha=0.3, marker='o',
                         norm=matplotlib.colors.PowerNorm(0.5),
                         s=(cat['area_exact'])**0.5)
    ax13.clear()
    pn = matplotlib.colors.PowerNorm(0.5, vmax=np.nanmax(cat['temperature_chi2']), vmin=np.nanmin(cat['temperature_chi2']))
    y = pn(cat['temperature_chi2'])
    colors = matplotlib.cm.jet(y)
    colors[:,3] = np.nan_to_num(np.array((pn(cat['temperature_chi2'])+1)/2.))
    sz = (cat['area_exact'])**0.5
    colors[:,3][sz>100] /= 3
    ln = matplotlib.colors.LogNorm(vmax=100, vmin=3, clip=True)
    colors[:,3] = np.nan_to_num(ln(sn))
    ax13.scatter(lon, cat['y_cen'],
                 c=colors,
                 marker='o',
                 edgecolor='none',
                 norm=matplotlib.colors.PowerNorm(0.5),
                 s=sz)
    ax13.set_xlabel("Galactic Longitude")
    ax13.set_ylabel("Galactic Latitude")
    cb13 = pl.colorbar(pts13)
    cb13.set_label(r'Temperature')


    try:
        brick = dendro.structure_at([262/(2 if 'smooth' in smooth else 1),143,725]).ancestor
        sgra = dendro.structure_at([239/(2 if 'smooth' in smooth else 1),97,907]).ancestor
        sgrb2 = dendro.structure_at([278/(2 if 'smooth' in smooth else 1),117,522]).ancestor
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
                   leaves_list=[sgra_leaves],
                   # r, b, g
                   color_list=['#CC4444', '#4444CC', '#44CC44'],
                   highlight_monotonic=True,
                   marker='s',
                   marker2=None,
                   linestyle='-', **kwargs):
        for leaves, color in zip(leaves_list,color_list):
            for leaf in leaves:
                xax,yax = [cat[leaf.idx][axname1]], [cat[leaf.idx][axname2]]
                axis.plot(xax, yax, marker, color=color, markeredgecolor='none', alpha=0.5)
                obj = leaf.parent
                while obj.parent:
                    xax.append(cat[obj.idx][axname1])
                    yax.append(cat[obj.idx][axname2])
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
    ax17.set_ylabel(r"Temperature [K]")
    fig17.savefig(fpath('dendrotem/temperature_vs_area{0}.pdf'.format(smooth)))

    fig17.clf()
    ax17 = fig17.gca()
    dendroplot(axis=ax17, axname1='Smean303', axname2='temperature_chi2', leaves_list=[sgra_leaves,sgrb2_leaves])
    ax17.set_xscale('log')
    ax17.set_xlabel(r"$\bar{S_\nu}(3_{03}-2_{02})$")
    ax17.set_ylabel(r"Temperature [K]")
    fig17.savefig(fpath('dendrotem/temperature_vs_flux{0}.pdf'.format(smooth)))

    fig18 = pl.figure(18)
    fig18.clf()
    ax18 = fig18.gca()
    dendroplot(leaves_list=[sgra_leaves, brick_leaves], axname2='temperature_chi2', axis=ax18)
    ax18.set_xscale('log')
    ax18.set_xlabel("Area (square arcseconds)")
    ax18.set_ylabel("Temperature [K]")
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

    hot = cat['temperature_chi2'] > 150
    ax23.errorbar(cat['v_rms'][hot], [149]*hot.sum(),
                  lolims=True, linestyle='none', capsize=0, alpha=alpha,
                  marker='^', color='r')
    for mask,color,alpha in masks_colors:
        ax23.errorbar(cat['v_rms'][mask], cat['temperature_chi2'][mask],
                    #yerr=[cat['elo_t'][mask], cat['ehi_t'][mask]],
                    linestyle='none', capsize=0, alpha=alpha, marker='.', color=color)
        ax23.set_xlabel("RMS Velocity")
        ax23.set_ylabel("Temperature [K]")
    ax23.set_ylim([0,150])

    # r,b,g = sgra, brick, sgrb2
    dendroplot(leaves_list=[sgra_leaves, brick_leaves, sgrb2_leaves],
               axname1='v_rms', axname2='temperature_chi2', axis=ax23,
               linestyle='none', marker='s', markerfacecolor='none',
               highlight_monotonic=False,)
    ax23.set_xlabel(r"RMS Velocity $(v_{rms}=\sigma_v)$")
    ax23.set_ylabel(r"Temperature [K]")
    fig23.savefig(fpath('dendrotem/dendro_temperature_vs_rmsvelocity{0}.pdf'.format(smooth)))



    for ii in range(1,13):
        pl.figure(ii)
        if ii not in (11,12,13,14,15):
            ax = pl.gca()
            ax.set_ylim(10, 125)
        pl.draw()

    pl.close(24)

    dview = dendro.viewer()
    structure = dendro.structure_at([262/(2 if 'smooth' in smooth else 1),143,725]).ancestor
    dview.hub.select(1, structure)
    dview.ax_image.axis((710,740,131,154))
    dview.slice_slider.set_val(262/(2 if 'smooth' in smooth else 1))
    dview.fig.savefig(fpath('dendrotem/dendrogram_viewer_brick{0}.pdf'.format(smooth)))

    

    pl.draw()
    pl.show()
