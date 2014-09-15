"""
Create an H2CO cloud catalog using dendrograms for further analysis using the temperature fitter
"""
from paths import hpath,mpath
from astrodendro import Dendrogram

def make_sncube(write=True):
    noise = fits.getdata(mpath('APEX_H2CO_merge_high_sub_noise.fits'))
    nhits = nhits = fits.getdata(mpath('APEX_H2CO_merge_high_nhits.fits'))
    noise[nhits<20] = np.nan

    sncube = fits.getdata(hpath('APEX_H2CO_303_202.fits')) / noise
    ff = fits.open(hpath('APEX_H2CO_303_202.fits'))
    ff[0].data = sncube

    if write:
        ff.writeto(hpath('APEX_H2CO_303_202_signal_to_noise_cube.fits'),
                   clobber=True)

    return sncube

def make_dend(sncube, view=True, write=True):

    dend = Dendrogram.compute(sncube, min_value=3, min_delta=2, min_npix=50,
                              verbose=True)

    if view:
        dend.viewer

    if write:
        dend.save_to(hpath("DendroMask_H2CO303202_signal_to_noise.hdf5"))

    return dend

if __name__ == "__main__":
    sncube = make_sncube()
    dend = make_dend(sncube)

