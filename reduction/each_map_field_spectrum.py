import numpy as np
import make_apex_cubes
import os
import pyspeckit
from pyspeckit.spectrum.readers import read_class
from astropy.time import Time


def make_spectra(dataset, maps):
    cl = read_class.ClassObject(os.path.join(make_apex_cubes.april2014path, dataset)+".apex")
    for map in maps:
        spectra = cl.get_spectra(sourcere='^'+map, telescope=list(cl.tels)[0])
        data = np.array([s for s,h in spectra])
        xarr = read_class.make_axis(header)
        sp = pyspeckit.Spectrum(data=data.mean(axis=0),
                                header=read_class.clean_header(header),
                                xarr=xarr)
        sp.write(paths.dpath('field_spectra/{name}_{window}_{date}.fits'.format(name=sp.specname,
                                                                                window=sp.header['XTEL'],
                                                                                date=Time(sp.header['DOBS'],
                                                                                          format='jyear').iso.split()[0])))


def make_all_spectra():
    for dataset in make_apex_cubes.datasets_2014:
        make_spectra(dataset, maps=make_apex_cubes.datasets_2014[dataset])
