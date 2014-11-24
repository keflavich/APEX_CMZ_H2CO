import numpy as np
import pylab as pl
import paths
import pyradex
import pyradex.fjdu
import warnings
warnings.filterwarnings('once')
from matplotlib.lines import Line2D
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))
pl.ion()


FF = pyradex.fjdu.Fjdu(species='ph2co-h2',
                   column=1e14,
                   temperature=50,
                   escapeProbGeom='lvg',
                   deltav=5,
                   collider_densities={'oH2':1e3, 'pH2':1e4})
RR = pyradex.Radex(species='ph2co-h2',
                   column=1e14,
                   temperature=50,
                   escapeProbGeom='lvg',
                   deltav=5,
                   collider_densities={'oH2':1e3, 'pH2':1e4})

drange1 = np.logspace(1,8,100)
fortho = 0.1

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()
for column, linestyle in zip((1e13,1e14,1e15),
                             (':','--','-')):
    for temperature,color in zip((20, 50, 100),
                                 ('r', 'b', 'purple')):
        tables1 = [RR(density={'oH2':density*fortho,
                               'pH2':density*(1-fortho)},
                      temperature=temperature,
                      column=column
                     )
                   for density in drange1]
        ax1.semilogx(drange1, [t['T_B'][12]/t['T_B'][2] for t in tables1],
                     label='{0} K'.format(temperature),
                     linestyle=linestyle,
                     color=color,
                    )
ax1.set_xlabel("Density ($n(\\mathrm{H}_2)$ cm$^{-3}$)")
ax1.legend((ax1.lines[-3], ax1.lines[-2], ax1.lines[-1],
            Line2D([0], [0], linestyle=':', color='k', linewidth=2), 
            Line2D([0], [0], linestyle='--', color='k', linewidth=2), 
            Line2D([0], [0], linestyle='-', color='k', linewidth=2), 
           ),
           ('$20$ K','$50$ K', '$100$ K',
            '$10^{13}$ cm$^{-2}$',
            '$10^{14}$ cm$^{-2}$',
            '$10^{15}$ cm$^{-2}$',),
            loc='best')
#ax1.set_title("$N(p-H_2CO)=10^{14}$")

fig1.savefig(paths.fpath('radex_ratio_lines_vs_density.pdf'), bbox_inches='tight')

crange1 = np.logspace(12,16.5,100)
fortho = 0.1

fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.gca()
for density, linestyle in zip((1e3,1e5,1e7),
                             (':','--','-')):
    for temperature,color in zip((20, 50, 100),
                                 ('r', 'b', 'purple')):
        tables1 = [RR(density={'oH2':density*fortho,
                               'pH2':density*(1-fortho)},
                      temperature=temperature,
                      column=column
                     )
                   for column in crange1]
        ax2.semilogx(crange1, [t['T_B'][12]/t['T_B'][2] for t in tables1],
                     label='{0} K'.format(temperature),
                     linestyle=linestyle,
                     color=color,
                    )
ax2.set_xlabel("Column ($N(\\mathrm{H}_2)$ cm$^{-2}$)")
ax2.legend((ax2.lines[-3], ax2.lines[-2], ax2.lines[-1],
            Line2D([0], [0], linestyle=':', color='k', linewidth=2), 
            Line2D([0], [0], linestyle='--', color='k', linewidth=2), 
            Line2D([0], [0], linestyle='-', color='k', linewidth=2), 
           ),
           ('$20$ K','$50$ K', '$100$ K',
            '$10^{3}$ cm$^{-3}$',
            '$10^{5}$ cm$^{-3}$',
            '$10^{7}$ cm$^{-3}$',),
            loc='best')
#ax2.set_title("$N(p-H_2CO)=10^{14}$")

fig2.savefig(paths.fpath('radex_ratio_lines_vs_column.pdf'), bbox_inches='tight')
