import numpy as np
import make_apex_cubes
import os
import paths
import pyspeckit
from pyspeckit.spectrum.readers import read_class
from astropy.time import Time


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

def fourier_flag_spectra(data, mode=3, fraction=0.05):
    """
    Fourier-based flagging: take the fft of the data along axis=1 (the spectral
    axis), then flag (ignore) the spectra with fourier `mode` in the top
    `fraction`
    """
    dft = np.fft.fft(data, axis=1)
    mode_data = np.abs(dft[:, mode])
    threshold = np.percentile(mode_data, 1-fraction)
    good = mode_data < threshold
    return data[good,:]

def make_all_spectra():
    for dataset in make_apex_cubes.datasets_2014:
        make_spectra(dataset, maps=make_apex_cubes.datasets_2014[dataset])
