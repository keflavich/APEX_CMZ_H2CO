import pylab as pl
import numpy as np
import aplpy
import os
import copy
from astropy import log
import paths
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))
from temperature_cubes import (tcube_dend, tcube_dend_smooth, tcube_direct,
                               tcubesm_direct)
from ratio_cubes import (ratiocube_303321, ratiocubesm_303321)
from masked_cubes import (cube303m, cube303msm)
from astropy import units as u
from image_registration.fft_tools import downsample

for ii in range(1,5):
    pl.close(ii)

if not os.path.isdir(paths.fpath('temvslon')):
    os.mkdir(paths.fpath('temvslon'))

figsize=(20,10)
vmin,vmax = -150.,150.
cbvmin,cbvmax = -80, 120
dv = 20.
vranges = np.arange(vmin, vmax, dv)
cm = pl.cm.rainbow_r

segmentdata = {'alpha': [(0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
               'blue': [(0.0, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
               'green': [(0.0, 0.0, 0.0), (0.5, 0.75, 0.75), (1.0, 0.0, 0.0)],
               'red': [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 1.0, 1.0)]}
cm = matplotlib.colors.LinearSegmentedColormap(name='rgb',
                                               segmentdata=segmentdata)

for weight in ("weighted_",""):
    for tcube,name in zip((tcube_dend, tcube_dend_smooth, tcube_direct,
                           tcubesm_direct, ratiocube_303321, ratiocubesm_303321),
                          ('dend','dendsm','direct','directsm','ratio','ratiosm')):

        log.info("Starting {0} {1}".format(name, weight))

        if weight:
            wcube = (cube303msm if 'sm' in name else cube303m)
            if tcube.shape != wcube.shape:
                log.info("Not weighting {0}".format(fn))
                continue
            weighted = copy.copy(tcube)
            weighted._data = wcube._data * tcube._data
            pv1 = weighted.sum(axis=1)
            pv2 = wcube.sum(axis=1)
            pv = pv1/pv2
            bmean = pv.value
        else:
            bmean = tcube.mean(axis=1)

        tcube_ds = tcube[::5 if 'sm' in name else 10,:,:] # for the WCS
        tcube_ds.data = downsample.downsample_axis(tcube._data,
                                                   5 if 'sm' in name else 10,
                                                   axis=0)

        bmeands = downsample.downsample_axis(bmean,
                                             5 if 'sm' in name else 10,
                                             axis=0)

        fig1 = pl.figure(1, figsize=figsize)
        fig1.clf()
        ax1 = fig1.gca()

        fig2 = pl.figure(2, figsize=figsize)
        fig2.clf()
        ax2 = fig2.gca()

        fig3 = pl.figure(3, figsize=figsize)
        fig3.clf()
        ax3 = fig3.gca()

        fig4 = pl.figure(4, figsize=figsize)
        fig4.clf()
        ax4 = fig4.gca()

        xax = tcube.world[0,0,:][2].value
        xax -= 360*(xax>180)
        specax = tcube_ds.spectral_axis.to(u.km/u.s).value

        assert specax.size == bmeands.shape[0]
        assert specax.size == tcube_ds.shape[0]

        gnoisescale = 5 if 'dend' in name else 0
        unoisescale = 1 if 'dend' in name else 0

        for v in vranges:
            color = cm((v-cbvmin)/(cbvmax-cbvmin))
            sel = (specax > v) & (specax < v+dv)
            x1 = np.argmin(np.abs(specax-v))
            x2 = np.argmin(np.abs(specax-(v+dv)))

            data = (tcube_ds.filled_data[x1:x2,:,:]
                    .T.reshape([xax.size,tcube_ds.filled_data[x1:x2,:,:].size/xax.size]))

            ax3.plot(xax,
                     data+np.random.randn(*data.shape)*gnoisescale,
                     marker=',', color=color,
                     linestyle='none', alpha=0.75)
            ax1.plot(xax, bmeands[sel,:].T, marker='.', color=color,
                     linestyle='none', alpha=0.75)
            ax2.plot(np.abs(xax), bmeands[sel,:].T, marker='.', color=color,
                     linestyle='none', alpha=0.75)
            ax4.plot(np.abs(xax),
                     data+(np.random.rand(*data.shape)*2-1)*unoisescale,
                     marker=',', color=color,
                     linestyle='none', alpha=0.75)

        for ax in (ax1,ax2,ax3,ax4):
            if 'ratio' in name:
                ax.set_ylim(0,0.60)
                ylabel = "Ratio $R_1$"
            else:
                ax.set_ylim(10,200)
                ylabel = "Temperature (K)"
        xlabel_lon = r"Galactic Longitude ($^\circ$)"
        xlabel_rad = r"$R_{projected}$ ($^{\circ}$)"
        cblabel = r"$v_{LSR}$ (km s$^{-1}$)"
        # http://stackoverflow.com/a/11558629/814354
        sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=cbvmin, vmax=cbvmax),
                                          cmap=cm)
        sm._A = []
        cb = fig1.colorbar(sm)
        ax1.set_xlim(1.6, -0.6)
        ax1.set_xlabel(xlabel_lon, labelpad=5)
        ax1.set_ylabel(ylabel, labelpad=5)
        cb.set_label(cblabel)

        ax3.set_xlim(1.6, -0.6)
        sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=cbvmin, vmax=cbvmax),
                                          cmap=cm)
        sm._A = []
        cb3 = fig3.colorbar(sm)
        ax3.set_xlabel(xlabel_lon, labelpad=5)
        ax3.set_ylabel(ylabel, labelpad=5)
        cb3.set_label(cblabel)

        ax2.set_xlim(0, 1.6)
        sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=cbvmin, vmax=cbvmax),
                                          cmap=cm)
        sm._A = []
        cb2 = fig2.colorbar(sm)
        ax2.set_xlabel(xlabel_rad, labelpad=5)
        ax2.set_ylabel(ylabel, labelpad=5)
        cb2.set_label(cblabel)


        ax4.set_xlim(0, 1.6)
        sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=cbvmin, vmax=cbvmax),
                                          cmap=cm)
        sm._A = []
        cb4 = fig4.colorbar(sm)
        ax4.set_xlabel(xlabel_rad, labelpad=5)
        ax4.set_ylabel(ylabel, labelpad=5)
        cb4.set_label(cblabel)

        log.info("Saving {0} {1}".format(name,weight))
        fig1.savefig(paths.fpath("temvslon/{1}temperature_vs_longitude_bmean_{0}.png".format(name,weight)))
        fig2.savefig(paths.fpath("temvslon/{1}temperature_vs_radius_bmean_{0}.png".format(name,weight)))
        fig3.savefig(paths.fpath("temvslon/{1}temperature_vs_longitude_{0}.png".format(name,weight)))
        fig4.savefig(paths.fpath("temvslon/{1}temperature_vs_radius_{0}.png".format(name,weight)))
        fig1.savefig(paths.fpath("temvslon/{1}temperature_vs_longitude_bmean_{0}.pdf".format(name,weight)))
        fig2.savefig(paths.fpath("temvslon/{1}temperature_vs_radius_bmean_{0}.pdf".format(name,weight)))
        fig3.savefig(paths.fpath("temvslon/{1}temperature_vs_longitude_{0}.pdf".format(name,weight)))
        fig4.savefig(paths.fpath("temvslon/{1}temperature_vs_radius_{0}.pdf".format(name,weight)))
        log.info("Completed {0} {1}".format(name,weight))
