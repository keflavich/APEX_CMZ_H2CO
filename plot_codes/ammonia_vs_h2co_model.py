from __future__ import print_function
from astropy import log
import numpy as np
import pyradex
import copy
import os
import pylab as pl
import paths
from astropy import table
from paths import analysispath
import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.utils.console import ProgressBar
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

density = 1e2
opr = 3
fortho = opr/(1.+opr)
temperatures = (50,100)

tex_nh311 = {}
tex_nh322 = {}
tex_h2co = {}
tex_h2co321 = {}

for temperature in temperatures:

    densities = np.logspace(1,9,100)
    tex_nh311[temperature] = []
    tex_nh322[temperature] = []
    tex_h2co[temperature] = []
    tex_h2co321[temperature] = []


    Rh2co = pyradex.Radex(species='ph2co-h2', column=1e12,
                          density={'oH2':density*fortho,
                                   'pH2':density*(1-fortho)},
                                    temperature=temperature)


    for density in ProgressBar(densities):

        Rh2co.density = {'oH2': density*fortho,
                         'pH2': density*(1-fortho),}
        Rh2co.temperature = temperature
        Rh2co.run_radex()

        tex_h2co[temperature].append(Rh2co.Tex[2])
        tex_h2co321[temperature].append(Rh2co.Tex[9])
        if np.any(Rh2co.tau[np.array([2,9])] > 0.5):
            print("Tau > 0.5 for T={0} n={1}".format(temperature,density))



    Rnh3 = pyradex.Radex(species='p-nh3', column=1e12,
                         density={'oH2':density*fortho, 'pH2':density*(1-fortho)},
                         temperature=temperature)


    print()
    for density in ProgressBar(densities):

        Rnh3.density = {'oH2': density*fortho,
                        'pH2': density*(1-fortho),}
        Rnh3.temperature = temperature
        Rnh3.run_radex()

        tex_nh311[temperature].append(Rnh3.Tex[8])
        tex_nh322[temperature].append(Rnh3.Tex[9])
        if np.any(Rnh3.tau[8:10] > 0.5):
            print("Tau > 0.5 for T={0} n={1}".format(temperature,density))
    print()

    tex_nh311[temperature] = np.array([x.value for x in tex_nh311[temperature]])
    tex_nh322[temperature] = np.array([x.value for x in tex_nh322[temperature]])
    tex_h2co[temperature] = np.array([x.value for x in tex_h2co[temperature]])
    tex_h2co321[temperature] = np.array([x.value for x in tex_h2co321[temperature]])

fig1 = pl.figure(1)
fig1.clf()
for axno in (1,2):
    ax = fig1.add_subplot(2,1,axno)
    linestyles = {100: '--', 20: '-', 50: ':'}
    for temperature in temperatures:
        ax.plot(densities, tex_nh311[temperature], color='b',
                linestyle=linestyles[temperature], alpha=0.5, linewidth=3,
                label="NH$_3$ 1-1 $T={0}$K".format(temperature))
        ax.plot(densities, tex_nh322[temperature], color='m',
                linestyle=linestyles[temperature], alpha=0.5, linewidth=3,
                label="NH$_3$ 2-2 $T={0}$K".format(temperature))
        ax.plot(densities, tex_h2co[temperature], color='r',
                linestyle=linestyles[temperature], alpha=0.5, linewidth=3,
                label="p-H$_2$CO $3_{{0,3}}-2_{{0,2}} T={0}$K".format(temperature))
        ax.plot(densities, tex_h2co321[temperature], color='g',
                linestyle=linestyles[temperature], alpha=0.5, linewidth=3,
                label="p-H$_2$CO $3_{{2,1}}-2_{{2,0}} T={0}$K".format(temperature))
    ax.vlines(1e3,0,55,color='k', alpha=0.5, linewidth=2, zorder=-10,
              linestyle='--')
    ax.vlines(1e4,0,55,color='k', alpha=0.5, linewidth=2, zorder=-10,
              linestyle=':')
    ax.vlines(1e5,0,55,color='k', alpha=0.5, linewidth=2, zorder=-10,
              linestyle='-')
    ax.set_xscale('log')
ax1,ax2 = fig1.axes
ax1.set_ylim(2.5,50)
ax2.set_ylim(2.5,5)
ax2.set_xlabel("Volume Density ($n(H_2)$ cm$^{-3}$)")
fig1.text(0.05, 0.5, "Excitation Temperature $T_{ex}$ (K)", rotation=90,
          ha='center', va='center')
fig1.subplots_adjust(hspace=0)
ax2.legend(loc='lower right', fontsize=14)
fig1.savefig(paths.fpath("radex/NH3andH2COexcitation.pdf"))

pl.draw()
pl.show()
