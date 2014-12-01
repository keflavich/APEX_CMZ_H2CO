import numpy as np
import pylab as pl
import paths
import pyradex
import pyradex.fjdu
import warnings
warnings.filterwarnings('once')
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))
pl.ion()


RR = pyradex.Radex(species='ph2co-h2',
                   column=1e14,
                   temperature=50,
                   escapeProbGeom='lvg',
                   deltav=5,
                   collider_densities={'oH2':1e3, 'pH2':1e4})

ff = 0.1
density = 1e5
fortho = 0.1
column = 1e14
temperatures = np.linspace(25,200,20)
temperatures = np.logspace(np.log10(25), np.log10(200), 40)
linewidth = 5.0

xarr = np.arange(-35,36,1, dtype='float')

noise = 0.1

tb303arr = []
tb321arr = []
recovratios = []
recovratiomeans = []
recovratiostds = []

for tem in temperatures:
    model_table = RR(collider_densities={'oH2':fortho*density,
                                         'pH2':(1-fortho)*density,},
                     column=column,
                     temperature=tem,
                     deltav=linewidth,
                    )
    tb303 = model_table['T_B'][2] *ff
    tb321 = model_table['T_B'][12]*ff

    line303 = tb303*np.exp(-xarr**2/(2*linewidth**2))
    line321 = tb321*np.exp(-xarr**2/(2*linewidth**2))
    int303 = line303.sum()
    int321 = line321.sum()

    recov303 = []
    recov321 = []
    recovratio = []
    for realization in np.arange(1500):
        noisearr = np.random.randn(xarr.size)
        line303n = line303 + noisearr*noise
        noisearr = np.random.randn(xarr.size)
        line321n = line321 + noisearr*noise

        selection = line321n > 2*noise
        if selection.sum() == 0:
            continue

        int303n = line303n[selection].sum()
        int321n = line321n[selection].sum()

        recovratio.append(int321n/int303n)

        recov303.append(int303n/int303)
        recov321.append(int321n/int321)

    recov_ratio_mean = np.mean(recovratio)
    recov_ratio_std = np.std(recovratio)

    tb303arr.append(tb303)
    tb321arr.append(tb321)
    recovratios.append(np.array(recovratio)/(tb321/tb303))
    recovratiomeans.append(recov_ratio_mean/(tb321/tb303))
    recovratiostds.append(recov_ratio_std)
            
from scipy.optimize import curve_fit
def f(x,p0,p1,p2):
    return p2-np.exp(-(x)**0.25*p1)*p0
p,c = curve_fit(f, tb321arr[2:], recovratiomeans[2:], p0=(10, 10, 1))

fig1=pl.figure(1)
fig1.clf()
ax1 = fig1.gca()
ax1.plot(tb321arr, f(np.array(tb321arr), *p), zorder=-5, alpha=0.5)
ax1.errorbar(tb321arr, recovratiomeans, yerr=recovratiostds, linestyle='none',
             marker='s')
ax1.set_xlabel("$T_B(3_{2,1}-2_{2,0})$", labelpad=20)
ax1.set_ylabel("Observed / Real $T(3_{2,1})/T(3_{0,3})$")
ax1.set_ylim(0.8, 2.0)
ax1.set_xlim(0.2, 1.6)
fig1.savefig(paths.fpath("noisemodeling_recoveredratio_means.pdf"), bbox_inches='tight')

fig2=pl.figure(2)
fig2.clf()
ax2 = fig2.gca()
ax2.boxplot(recovratios, positions=tb321arr, widths=0.010)
ax2.set_xticks(tb321arr[::4])
ax2.set_xticklabels(["{0:0.2f}".format(x) for x in tb321arr[::4]])
ax2.set_xlabel("$T_B(3_{2,1}-2_{2,0})$", labelpad=20)
ax2.set_ylabel("Observed / Real $T(3_{2,1})/T(3_{0,3})$")
ax2.grid()
ax2.set_ylim(0.8, 2.0)
ax2.set_xlim(0.2, 1.6)
fig2.savefig(paths.fpath("noisemodeling_recoveredratio_boxplots.pdf"), bbox_inches='tight')

# Not useful: had an idea but it didn't work out.  Save for later...
#fig3=pl.figure(2)
#fig3.clf()
#ax3 = fig3.gca()
#ax3.boxplot(recovratios, positions=recovratiomeans, widths=0.010)
#ax3.set_xticks(tb321arr[::4])
#ax3.set_xticklabels(["{0:0.2f}".format(x) for x in tb321arr[::4]])
#ax3.set_xlabel("$T_B(3_{2,1}-2_{2,0})$", labelpad=20)
#ax3.set_ylabel("Observed / Real $T(3_{2,1})/T(3_{0,3})$")
#ax3.grid()
#ax3.set_ylim(0.8, 2.0)
#ax3.set_xlim(0.2, 1.6)
#fig3.savefig(paths.fpath("noisemodeling_recoveredratio_boxplots.pdf"), bbox_inches='tight')
