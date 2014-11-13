import numpy as np
import os
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel

def suppress_frange(data, frange, convinterp=True, width=2, convolve_with_fft=False):
    """
    Given a spectrum, interpolate across a specified frequency range.
    Frequency is in 1/pixel units.

    Parameters
    ----------
    data : np.ndarray
        One dimensionala array containing the spectrum to be squished
    frange : (min, max)
        The range of frequencies (in 1/pixel units) to suppress
    convinterp : bool
        Use convolution for the interpolation?  Otherwise, does simple linear
        interpolation (which can be unstable)
    width : int
        Width of the convolution kernel
    convolve_with_fft : bool
        Use `astropy.convolution.convolve` or `astropy.convolution.convolve_fft`.
        For 1D arrays like this, convolve is generally faster.
    """

    ft = np.fft.fft(data)
    freq = np.fft.fftfreq(data.size)

    lmask = (abs(freq) < frange[1])
    umask = (abs(freq) > frange[0])
    mask = lmask & umask

    # Select one side of the data
    midpt = len(mask)/2
    whmask, = np.where(mask)
    startpt = whmask[whmask<midpt][0] - 1
    endpt = whmask[whmask<midpt][-1] + 1

    if endpt >= len(mask):
        endpt = len(mask)-1

    if convinterp:
        mdata = ft.copy()
        mdata[mask] = np.nan
        kernel = Gaussian1DKernel(width, x_size=width*8+1)
        kernel.normalize()
        if convolve_with_fft:
            amp = convolve_fft(np.abs(mdata), kernel, interpolate_nan=True)
            phase = convolve_fft(np.angle(mdata), kernel, interpolate_nan=True)
        else:
            amp = convolve(np.abs(mdata), kernel, boundary='extend')
            phase = convolve(np.angle(mdata), kernel, boundary='extend')
        interpdata = (np.cos(phase)*amp+1j*np.sin(phase)*amp)[mask]
        ft[mask] = interpdata
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

    return np.fft.ifft(ft).real
