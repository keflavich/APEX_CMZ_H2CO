from individual_spectra import *
import full_cubes
import paths

spd = full_cubes.cube_merge_high.mean(axis=(1,2))
sp = pyspeckit.Spectrum(xarr=full_cubes.cube_merge_high.spectral_axis, data=spd)
sp.xarr.refX = full_cubes.pcube_merge_high.xarr.refX
sp.xarr.refX_units = full_cubes.pcube_merge_high.xarr.refX_units
sp.specname = 'Whole CMZ'
sp.unit = 'K'

sp.plotter()
sp.specfit.Registry.add_fitter('h2co_simple', simple_fitter2, 6,
                               multisingle='multi')
sp.baseline(exclude=[151, 761, 916, 1265], selectregion=True,
            highlight_fitregion=True, xtype='pixel', subtract=True)

sp.error[:] = sp.data[761:916].std()

sp.plotter()
sp.specfit(fittype='h2co_simple', multifit=True,
           guesses=[0.06, 10, 20, 0.5, 0.7, 0.03],
           limited=[(True,True)] * 6,
           limits=[(0,20),[-150,150],(1, 60),(0,1),(0.3,1.1),(0,1e5)],
          )

sp.specfit.plot_fit(show_components=True)
sp.plotter.savefig(paths.fpath('simple/WholeCMZ_6parameter.pdf'),
                   bbox_inches='tight')
