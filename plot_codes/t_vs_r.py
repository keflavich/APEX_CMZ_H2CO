import pylab as pl
import numpy as np
import aplpy
import os
import copy
from astropy import log
import paths
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))
from temperature_cubes import tcube_dend, tcube_dend_smooth, tcube_direct

if not os.path.isdir(paths.fpath('temvslon')):
    os.mkdir(paths.fpath('temvslon'))

vmin,vmax = -150.,190.
dv = 20.
vranges = np.arange(vmin, vmax, dv)
cm = pl.cm.rainbow_r

segmentdata = {'alpha': [(0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
               'blue': [(0.0, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
               'green': [(0.0, 0.0, 0.0), (0.5, 0.75, 0.75), (1.0, 0.0, 0.0)],
               'red': [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 1.0, 1.0)]}
cm = matplotlib.colors.LinearSegmentedColormap(name='rgb',
                                               segmentdata=segmentdata)

for tcube,name in zip((tcube_dend, tcube_dend_smooth), # , tcube_direct
                      ('dend','dendsm')):#,'direct'

    log.info("Starting {0}".format(name))

    bmean = tcube.median(axis=1)

    fig1 = pl.figure(1)
    fig1.clf()
    ax1 = fig1.gca()

    xax = tcube.world[0,0,:][2].value
    xax -= 360*(xax>180)
    specax = tcube.spectral_axis.to(u.km/u.s).value

    for v in vranges:
        color = cm((v-vmin)/(vmax-vmin))
        sel = (specax > v) & (specax < v+dv)
        x1 = np.argmin(np.abs(specax-v))
        x2 = np.argmin(np.abs(specax-(v+dv)))

        #ax1.plot(xax,
        #         (tcube.filled_data[x1:x2,:,:]
        #          .T.reshape([xax.size,tcube.filled_data[x1:x2,:,:].size/xax.size])),
        #         marker=',', color=color,
        #         linestyle='none', alpha=0.75)
        ax1.plot(xax, bmean[sel,:].T, marker=',', color=color,
                 linestyle='none', alpha=0.75)

    ax1.set_xlim(1.6, -0.6)
    ax1.set_ylim(20,150)
    # http://stackoverflow.com/a/11558629/814354
    sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax),
                                      cmap=cm)
    sm._A = []
    cb = fig1.colorbar(sm)
    ax1.set_xlabel(r"Galactic Longitude ($^\circ$)", labelpad=5)
    ax1.set_ylabel("Temperature (K)", labelpad=5)
    cb.set_label(r"$v_{LSR}$ (km s$^{-1}$)")

    log.info("Saving {0}".format(name))
    fig1.savefig(paths.fpath("temvslon/temperature_vs_longitude_{0}.pdf".format(name)))
    log.info("Completed {0}".format(name))
