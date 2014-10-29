"""
Perform fourier analysis (i.e., look for frequency peaks) on the PCA-extracted
Most Correlated Components
"""
from make_apex_cubes import june2013datapath, april2014path
import numpy as np
from astropy.io import fits
import os
import pylab as pl
import paths
import matplotlib
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
import pyspeckit
matplotlib.rc_file(paths.pcpath('pubfiguresrc'))

def find_pca_files():
    pass


def spline_removal(fn, linecen=1500, linewidth=20, lineamp=10):
    sp = pyspeckit.Spectrum(fn)
    sp.plotter()
    sp.baseline(spline=True, subtract=False, spline_sampling=500, order=3)
    bl1 = sp.baseline.basespec
    synthline = lineamp*np.exp(-(sp.xarr-linecen)**2/(2*linewidth**2.))
    sp2 = sp.copy()
    sp2.data += synthline
    sp2.plotter()
    sp2.baseline.set_spectofit()
    sp2.baseline(spline=True, subtract=False, spline_sampling=500, order=3)
    return sp,sp2

def do_example_fsupp(interpolate=False, vrange=[100,300]):
    """
    Run a specific example on a specific extracted PCA component
    """
    fn = os.path.join(june2013datapath,
                      'M-091.F-0019-2013-2013-06-12/AP-H201-X202_pca_component_0.fits')
    fourier_suppression(fn, vrange=vrange, save=True, linecen=-1500.,
                        interpolate=interpolate)
    return fn

def fourier_suppression(fn, withline=False, vrange=[100,'max'], linewidth=20.,
                        linecen=250., save=False, suppression_factor=3,
                        interpolate=False, convinterp=True):
    """
    Given a spectrum (fn), suppress a velocity range `vrange` by dividing it by
    some factor.  Then plot...
    """

    e1 = fits.getdata(fn)

    hdr = fits.getheader(fn)
    vres = hdr['CDELT1']
    ff = np.fft.fftfreq(e1.size)

    fvr = e1.size*vres
    fvelo = fvr/(ff/ff[1])/2.
    fvelo[0] = fvr
    velo = vres * (np.arange(e1.size)+1-hdr['CRPIX1'])+hdr['CRVAL1']

    # Add a 20 km/s wide Gaussian line, see what happens to it
    line = np.exp(-(velo+linecen)**2/(linewidth**2*2.)) * 10
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


    if interpolate:
        midpt = len(mask)/2
        whmask, = np.where(mask)
        startpt = whmask[whmask<midpt][0] - 1
        endpt = whmask[whmask<midpt][-1] + 1
        if endpt >= len(mask):
            endpt = len(mask)-1
        for fff in (ft,ft2):

            if convinterp:
                mdata = fff.copy(); mdata[mask] = np.nan
                #realdata = convolve(mdata.real, Gaussian1DKernel(1, x_size=51), boundary='extend')
                #imagdata = convolve(mdata.imag, Gaussian1DKernel(1, x_size=51), boundary='extend')
                #interpdata = (realdata[mask]+1j*imagdata[mask]) 
                amp = convolve(np.abs(mdata), Gaussian1DKernel(1, x_size=51), boundary='extend')
                phase = convolve(np.angle(mdata), Gaussian1DKernel(1, x_size=51), boundary='extend')
                interpdata = (np.cos(phase)*amp+1j*np.sin(phase)*amp)[mask]
                #pl.figure(5)
                #pl.clf()
                #ax1 = pl.subplot(3,1,1)
                #ax1.plot(np.arange(1,30), fff.real[1:30])
                #ax1.plot(np.arange(1,30), mdata.real[1:30])
                #ax1.plot(whmask[:len(whmask)/2], interpdata.real[:len(whmask)/2],'--')
                #ax2 = pl.subplot(3,1,2)
                #ax2.plot(np.arange(1,30), fff.imag[1:30])
                #ax2.plot(np.arange(1,30), mdata.imag[1:30])
                #ax2.plot(whmask[:len(whmask)/2], interpdata.imag[:len(whmask)/2],'--')
                #ax3 = pl.subplot(3,1,3)
                #ax3.plot(np.arange(1,30), abs(fff)[1:30])
                #ax3.plot(np.arange(1,30), abs(mdata)[1:30])
                #ax3.plot(whmask[:len(whmask)/2], abs(interpdata[:len(whmask)/2]),'--')
                fff[mask] = interpdata
                #ax1.plot(np.arange(1,30), fff.real[1:30],':')
                #ax2.plot(np.arange(1,30), fff.imag[1:30],':')
                #ax3.plot(np.arange(1,30), abs(fff)[1:30],':')
            else:
                for order,compare in zip((-1,1),(np.less,np.greater_equal)):
                    mm = mask & compare(np.arange(len(mask), dtype='int'), midpt)
                    realdata = np.interp(np.arange(startpt+1, endpt),
                                              [startpt,endpt],
                                              [fff.real[startpt],fff.real[endpt]])
                    imagdata = np.interp(np.arange(startpt+1, endpt),
                                              [startpt,endpt],
                                              [fff.imag[startpt],fff.imag[endpt]])
                    fff[mm] = (realdata+1j*imagdata)[::order]
    else:
        level = max(abs(ft)[levelmask])/suppression_factor
        level2 = max(abs(ft2)[levelmask])/suppression_factor

        ft[(mask)] /= abs(ft[(mask)])/level
        ft2[(mask)] /= abs(ft2[(mask)])/level2

    ax1 = pl.figure(1).gca()
    ax1.cla()
    ax1.plot(velo, e1, linestyle='none', color='k', marker=',', label='Input')
    ax1.plot(velo, np.fft.ifft(ft)+offset, alpha=0.8, linestyle='none',
             marker=',', color='r', label='FT suppressed')
    ax1.set_xlabel("Velocity")
    ax1.set_ylabel("Brightness (K-ish)")
    leg = ax1.legend(loc='best', markerscale=20)
    for ll in leg.get_lines():
        ll.set_marker('o')
        ll.set_markersize(5)
    if save:
        pl.figure(1).savefig(paths.fpath('baselines/ft_suppression_spectra.png'))

    ax4 = pl.figure(4).gca()
    ax4.cla()
    ax4.set_title("Synthetic {0} km/s line".format(linewidth))
    ax4.plot(velo, e2-np.median(e2), linestyle='none', color='b', marker=',', zorder=-1,
             alpha=0.5, label='Input')
    ax4.plot(velo, np.fft.ifft(ft2)+offset-np.median(e2), alpha=0.5, linestyle='none',
             marker=',', color='g', zorder=-1, label='FT Suppressed')
    ax4.plot(velo, e2-np.fft.ifft(ft2)-offset, alpha=0.5, linestyle='none',
             marker=',', color='r', zorder=-1, label='(a) Diff')
    ax4.plot(velo, e1-np.fft.ifft(ft)-offset*2, alpha=0.5, linestyle='none',
             marker=',', color='m', zorder=-1, label='(b) Diff with no synthetic line')
    ax4.plot(velo, (e2-np.fft.ifft(ft2))-(e1-np.fft.ifft(ft))-offset*3, alpha=0.5, linestyle='none',
             marker=',', color='k', zorder=-1, label='(a)-(b)')
    ax4.plot(velo, line+offset, alpha=0.5,
             linewidth=0.5, color='k', zorder=1, label='Synthetic Line')
    ax4.set_xlabel("Velocity")
    ax4.set_ylabel("Brightness (K-ish)")
    leg = ax4.legend(loc='lower left', markerscale=20, fontsize=16)
    for ll in leg.get_lines():
        ll.set_marker('o')
        ll.set_markersize(5)
    if save:
        pl.figure(4).savefig(paths.fpath('baselines/ft_suppression_spectra_synthline.png'))
        ax4.axis([1300,1700,-35,21])
        pl.figure(4).savefig(paths.fpath('baselines/ft_suppression_spectra_synthline_zoom.png'))

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
    leg = ax3.legend(loc='best', markerscale=20)
    if save:
        pl.figure(3).savefig(paths.fpath('baselines/ft_suppression_differences.png'))


    pl.draw()
    pl.show()
