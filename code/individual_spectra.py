import pyspeckit
from paths import h2copath, mergepath, figurepath
import os
from pyspeckit_fitting import h2co_radex_fitter, simplemodel, simple_fitter

cube = pyspeckit.Cube(os.path.join(mergepath, 'APEX_H2CO_merge_high.fits'))

spectra = {
        'brick1': cube.get_apspec((0.2426,0.0081,30), coordsys='galactic', wunit='arcsec'),
        'brick2': cube.get_apspec((0.2583,0.0181,30), coordsys='galactic', wunit='arcsec'),
    }

pars = {
    'brick1': {'ncomp': 1},
    'brick2': {'ncomp': 2},
}

for spname in ('brick1','brick2'):
    ncomp = pars[spname]['ncomp']
    sp = spectra[spname]
    sp.specname = spname
    sp.specfit.Registry.add_fitter('h2co_mm_radex', h2co_radex_fitter, 5,
                             multisingle='multi')
    sp.specfit.Registry.add_fitter('h2co_simple', simple_fitter, 4, multisingle='multi')
    sp.error[:] = sp.stats([218e9, 218.1e9])['std']
    sp.specfit(fittype='h2co_simple', multifit=True,
               guesses=[1,25,5,0.5,1]*ncomp)

    sp.plotter()
    sp.specfit.plot_fit()
    sp.plotter.savefig(os.path.join(figurepath,
                                    "{0}_fit_4_lines_simple.pdf".format(spname)))

    sp.specfit(fittype='h2co_mm_radex', multifit=True,
               guesses=[100,14,4.5,
                        sp.specfit.parinfo['VELOCITY0'].value,
                        sp.specfit.parinfo['WIDTH0'].value]*ncomp,
               limits=[(20,200),(11,15),(3,5.5),(-105,105),(1,8)]*ncomp,
               limited=[(True,True)]*5*ncomp,
               fixed=[False,False,False,True,True]*ncomp,
               quiet=False,)
    sp.plotter.savefig(os.path.join(figurepath,
                                    "{0}_fit_h2co_mm_radex.pdf".format(spname)))
