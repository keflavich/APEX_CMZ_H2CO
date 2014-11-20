import numpy as np
import make_apex_cubes
import os
import paths
import pyspeckit
from pyspeckit.spectrum.readers import read_class
from astropy.time import Time
import pylab as pl


def make_spectra(dataset, maps, tels=['AP-H201-X202', 'AP-H201-X201']):
    cl = read_class.ClassObject(os.path.join(make_apex_cubes.april2014path, dataset)+".apex")
    for tel in tels:
        for map in maps:
            spectra = cl.get_spectra(sourcere='^'+map, telescope=tel)
            data = np.array([s for s,h in spectra])
            header = spectra[0][1]
            xarr = read_class.make_axis(header)
            sp = pyspeckit.Spectrum(data=data.mean(axis=0),
                                    header=read_class.clean_header(header),
                                    xarr=xarr)
            sp.write(paths.dpath('field_spectra/{name}_{window}_{date}.fits'.format(name=sp.specname,
                                                                                    window=sp.header['XTEL'],
                                                                                    date=Time(sp.header['DOBS'],
                                                                                              format='jyear').iso.split()[0])))
            
            datafflag = fourier_flag_spectra(data)

            header['FTFLGMOD'] = (3, 'Fourier component flagged')
            header['FTFLGFRC'] = (0.05, 'Bad fraction of Fourier component')
            sp = pyspeckit.Spectrum(data=datafflag.mean(axis=0),
                                    header=read_class.clean_header(header),
                                    xarr=xarr)
            sp.write(paths.dpath('field_spectra/{name}_{window}_{date}_fflagged.fits'.format(name=sp.specname,
                                                                                             window=sp.header['XTEL'],
                                                                                             date=Time(sp.header['DOBS'],
                                                                                                       format='jyear').iso.split()[0])))

def make_fourier_spectra(dataset, maps, tels=['AP-H201-X202', 'AP-H201-X201']):
    for tel in tels:
        for map in maps:
            sp = pyspeckit.Spectrum(paths.dpath('field_spectra/'
                                                "{name}_{window}_{date}{ff}.fits".format(name=map,
                                                                                         window=tel,
                                                                                         date=dataset[-10:],
                                                                                         ff="")))
            spff = pyspeckit.Spectrum(paths.dpath('field_spectra/'
                                                  "{name}_{window}_{date}{ff}.fits".format(name=map,
                                                                                          window=tel,
                                                                                          date=dataset[-10:],
                                                                                          ff="_fflagged")))

            ft1 = np.fft.fft(sp.data)
            ft2 = np.fft.fft(spff.data)

            ftf = np.fft.fftfreq(sp.data.size)
            vres = float(sp.header['VRES'])
            fvr = ftf.size*vres
            fvelo = fvr/(ftf/ftf[1])/2.

            fig = pl.figure(1)
            fig.clf()
            ax = fig.gca()
            ax2 = ax.twiny()
            def tick_function(X):
                V = abs(fvr/(X/ftf[1])/2.)
                return ["%0.3g" % z for z in V]

            pl.loglog(ftf, np.abs(ft1), 'k', alpha=0.5)
            pl.loglog(ftf, np.abs(ft2), 'r', alpha=0.5)

            tick_locs = np.logspace(-np.log10(ftf.size/2.), np.log10(0.5), 7)
            ax2.set_xticks(tick_locs)
            ax2.set_xticklabels(tick_function(tick_locs))
            ax.set_xlabel("Fourier Frequency")
            ax2.set_xlabel("Velocity Scale (km/s)")
            ax.set_ylabel("Fourier Amplitude")

            ax.set_xlim(ftf[1], ftf[ftf.size/2-1])
            ax.set_xscale('log')
            ax.set_xlim(ftf[1], ftf[ftf.size/2-1])
            ax.set_xlim(ftf[1], ftf[ftf.size/2-1])
            ax.set_xlim(ftf[1], ftf[ftf.size/2-1])
            pl.xlim(ftf[1], ftf[ftf.size/2-1])
            assert ax.get_xlim() == (ftf[1], ftf[ftf.size/2-1])

            pl.savefig(paths.dpath('field_spectra/'
                                   "{name}_{window}_{date}_psds.png".format(name=map,
                                                                            window=tel,
                                                                            date=dataset[-10:],)
                                  ),
                       bbox_inches='tight')


def fourier_flag_spectra(data, mode=3, fraction=0.05):
    """
    Fourier-based flagging: take the fft of the data along axis=1 (the spectral
    axis), then flag (ignore) the spectra with fourier `mode` in the top
    `fraction`
    """
    dft = np.fft.fft(data, axis=1)
    mode_data = np.abs(dft[:, mode])
    threshold = np.percentile(mode_data, 100*(1-fraction))
    good = mode_data < threshold
    return data[good,:]

def make_all_spectra():
    for dataset in make_apex_cubes.datasets_2014.keys()[::-1]:
        make_spectra(dataset, maps=make_apex_cubes.datasets_2014[dataset])

def make_all_fourier_spectra():
    pl.ioff()
    for dataset in make_apex_cubes.datasets_2014.keys()[::-1]:
        make_fourier_spectra(dataset, maps=make_apex_cubes.datasets_2014[dataset])
    pl.ion()
