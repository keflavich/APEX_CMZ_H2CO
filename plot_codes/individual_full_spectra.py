import pyregion
import numpy as np
import pyspeckit
from astropy import table
import spectral_cube
from spectral_cube import SpectralCube
from paths import (h2copath, mergepath, figurepath, regpath, analysispath,
                   mpath, dpath)
import os
from astropy.io import fits
import photutils
from astropy import coordinates
from astropy import units as u
from astropy import log
from astropy.utils.console import ProgressBar
import pylab as pl
from astroquery.splatalogue import Splatalogue

from shfi_otf_pipeline.make_apex_cubes import all_lines

pl.ioff()
pl.figure(1, figsize=(10,6))

def fpath(x, figurepath=os.path.join(figurepath, 'fullspectra')):
    return os.path.join(figurepath, x)

regs = pyregion.open(regpath+'spectral_apertures.reg')

ftemplate =  'APEX_H2CO_2014_merge_{0}.fits'

# Nice idea, but Splatalogue doesn't know what lines are really there
#line_table = Splatalogue.query_lines(216.9*u.GHz,
#                                     221*u.GHz, top20='top20',
#                                     energy_max=150, energy_type='eu_k',
#                                     line_lists=['SLAIM'])
#names = np.array(['{0}_{1}'.format(a,b) for a,b in zip(line_table['Species'],
#                                                       line_table['Resolved QNs'])])
#lfreq = line_table['Freq-GHz']

for lh in ('low','high'):
    cube = SpectralCube.read(mpath(ftemplate.format(lh)))
    noisehdu = cube.std(axis=0).hdu
    noisehdr = noisehdu.header
    noise = noisehdu.data
    noiseokmask = np.isfinite(noise)

    xarr = pyspeckit.units.SpectroscopicAxis(cube.spectral_axis.value,
                                             unit=str(cube.spectral_axis.unit),
                                             refX=cube.wcs.wcs.restfrq,
                                             refX_units='Hz')

    spectra = {}
    for region_number,reg in enumerate(regs):
        name = reg.attr[1]['text']
        if name not in spectra:
            #sp = cube.get_apspec(reg.coord_list,coordsys='galactic',wunit='degree')
            shape = pyregion.ShapeList([reg])
            mask = shape.get_mask(header=noisehdr, shape=noise.shape)
            scube = cube.subcube_from_ds9region(shape)
            data = scube.apply_numpy_function(np.nanmean, axis=(1,2))
            error = ((noise[mask & noiseokmask]**2).sum()**0.5/np.count_nonzero(mask))
            sp = pyspeckit.Spectrum(data=data,
                                    error=np.ones(data.size)*error,
                                    xarr=xarr, header=cube.wcs.to_header())
            sp.xarr.convert_to_unit('GHz')
            sp.header['ERROR'] = error
            #sp.error[:] = sp.stats((218.5e9,218.65e9))['std']
            sp.specname = reg.attr[1]['text']
            # Error is already computed above; this is an old hack
            #sp.error[:] = sp.stats((218e9,218.1e9))['std']
            spectra[name] = sp
            sp.unit = "$T_{A}$ [K]"
        else:
            sp = spectra[name]

        sp.plotter.figure = pl.figure(1)
        sp.plotter(errstyle='fill')
        #linesel = (lfreq > sp.xarr.as_unit('GHz').min()) & (lfreq < sp.xarr.as_unit('GHz').max())
        try:
            #sp.plotter.line_ids(names[linesel], lfreq[linesel], xval_units='GHz')
            sp.plotter.line_ids(all_lines.keys(), all_lines.values(), xval_units='GHz')
        except RuntimeError as ex:
            print ex

        spname = sp.specname.replace(" ","_")
        sp.plotter.savefig(fpath("{0}_{1}.png".format(spname, lh)))
