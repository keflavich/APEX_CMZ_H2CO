"""
Perform fourier analysis (i.e., look for frequency peaks) on the PCA-extracted
Most Correlated Components
"""
from make_apex_cubes import june2013datapath, april2014path
import numpy as np
from astropy.io import fits
import pylab as pl
import paths
import matplotlib
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

def find_pca_files():
    pass

def something():
    for sp in spectra:
        ft = np.fft.fft(sp.data)
        ff = np.fft.fftfreq(sp.data.size)
        pos = ff >= 0
        fig = pl.figure(1)
        ax = fig.gca()
        ax.loglog(ff[pos], abs(ft[pos]))

def fourier_suppression(fn, withline=False, vrange=[100,'max'], linewidth=20.,
                        save=False):

    e1 = fits.getdata(fn)

    hdr = fits.getheader(fn)
    vres = hdr['CDELT1']
    ff = np.fft.fftfreq(e1.size)

    fvr = e1.size*vres
    fvelo = fvr/(ff/ff[1])/2.
    fvelo[0] = fvr
    velo = vres * (np.arange(e1.size)+1-hdr['CRPIX1'])+hdr['CRVAL1']

    # Add a 20 km/s wide Gaussian line, see what happens to it
    line = np.exp(-(velo+250.)**2/(linewidth**2*2.)) * 10
    e2 = line+e1

    ft = np.fft.fft(e1)
    ft2 = np.fft.fft(e2)


    #ft[(abs(fvelo)>100) & (abs(fvelo<165))] /= abs(ft[(abs(fvelo)>100) & (abs(fvelo<165))])/1000
    #ft[(abs(fvelo)>230) & (abs(fvelo<556))] /= abs(ft[(abs(fvelo)>230) & (abs(fvelo<556))])/1000
    if 'max' in vrange:
        vrange[1] = max(fvelo)+1
        offset = 0
    else:
        offset = 10

    lmask = (abs(fvelo) < vrange[1])
    umask = (abs(fvelo) > vrange[0])
    mask = lmask & umask
    levelmask = (abs(fvelo) < vrange[0])

    level = max(abs(ft)[levelmask])/3

    level2 = max(abs(ft2)[levelmask])/3

    ft[(mask)] /= abs(ft[(mask)])/level
    ft2[(mask)] /= abs(ft2[(mask)])/level2

    ax1 = pl.figure(1).gca()
    ax1.cla()
    ax1.plot(velo, e1, linestyle='none', color='k', marker=',', label='Input')
    ax1.plot(velo, np.fft.ifft(ft)+offset, alpha=0.8, linestyle='none',
             marker=',', color='r', label='FT suppressed')
    ax1.set_xlabel("Velocity")
    ax1.set_ylabel("Brightness (K-ish)")
    ax1.legend(loc='best')
    if save:
        pl.figure(1).savefig(paths.fpath('baselines/ft_suppression_spectra.png'))

    ax4 = pl.figure(4).gca()
    ax4.cla()
    ax4.set_title("Synthetic {0} km/s line".format(linewidth))
    ax4.plot(velo, e2, linestyle='none', color='b', marker=',', zorder=-1,
             alpha=0.5, label='Input')
    ax4.plot(velo, np.fft.ifft(ft2)+offset, alpha=0.5, linestyle='none',
             marker=',', color='g', zorder=-1, label='FT Suppressed')
    ax4.set_xlabel("Velocity")
    ax4.set_ylabel("Brightness (K-ish)")
    ax4.legend(loc='best')
    if save:
        pl.figure(4).savefig(paths.fpath('baselines/ft_suppression_spectra_synthline.png'))

    ax2 = pl.figure(2).gca()
    ax2.cla()
    ax2.loglog(fvelo, abs(np.fft.fft(e1)), linewidth=1, label='Input')
    ax2.loglog(fvelo, abs(np.fft.fft(e2)), linewidth=1,
               label='Synthetic {0} km/s line'.format(linewidth))
    ax2.loglog(fvelo, abs(ft), linewidth=1, label='FT Suppressed (Input)')
    ax2.loglog(fvelo, abs(ft2), linewidth=1, label='FT Suppressed (Synthetic)')
    ax2.set_xlabel("Velocity Scale")
    ax2.set_ylabel("$|FT(spectrum)|$")
    ax2.legend(loc='best')
    if save:
        pl.figure(2).savefig(paths.fpath('baselines/ft_suppression_fourierpower.png'))

    ax3 = pl.figure(3).gca()
    ax3.cla()
    ax3.plot(velo, e1-np.fft.ifft(ft), label='Original')
    ax3.plot(velo, e2-np.fft.ifft(ft2), label='Synthetic {0} km/s line'.format(linewidth))
    ax3.set_xlabel("Velocity")
    ax3.set_ylabel("Difference (Input-FT_suppressed)")
    ax3.legend(loc='best')
    if save:
        pl.figure(3).savefig(paths.fpath('baselines/ft_suppression_differences.png'))


    pl.draw()
    pl.show()
