#from individual_spectra import *
import numpy as np
import pyspeckit
from astropy import units as u
from pyspeckit_fitting import simple_fitter2
import full_cubes
import paths
import pylab as pl

full_cubes.cube_merge_high.allow_huge_operations = True
spd = full_cubes.cube_merge_high.mean(axis=(1,2)).value
sp = pyspeckit.Spectrum(xarr=full_cubes.cube_merge_high.spectral_axis,
                        data=spd,
                        xarrkwargs={'velocity_convention': 'radio'},
                       )
sp.xarr.refX = full_cubes.pcube_merge_high.xarr.refX
sp.xarr.refX_unit = full_cubes.pcube_merge_high.xarr.refX_unit
sp.xarr.convert_to_unit(u.GHz)
sp.specname = 'Whole CMZ'
sp.unit = 'K'

sp.plotter()
sp.specfit.Registry.add_fitter('h2co_simple', simple_fitter2, 6,
                               multisingle='multi')
sp.baseline(exclude=[151, 761, 916, 1265], selectregion=True,
            highlight_fitregion=True, xtype='pixel', subtract=True)

sp.error[:] = sp.data[761:916].std()

sp.plotter()
sp.specfit(fittype='h2co_simple',
           guesses=[0.06, 10, 20, 0.5, 0.7, 0.03],
           limited=[(True,True)] * 6,
           limits=[(0,20),[-150,150],(1, 60),(0,1),(0.3,1.1),(0,1e5)],
          )

sp.specfit.plot_fit(show_components=True)
sp.plotter.savefig(paths.fpath('simple/WholeCMZ_6parameter.pdf'),
                   bbox_inches='tight')

lat,lon = full_cubes.cube_merge_high.world[0,:,:][1:]
sgrb2_cloud_mask = ((lon-0.674*u.deg)**2 + (lat+0.027*u.deg)**2)**0.5 < 5*u.arcmin
nob2cube = full_cubes.cube_merge_high.with_mask(~sgrb2_cloud_mask)
spd_nob2 = nob2cube.mean(axis=(1,2)).value

print("Sgr B2 pix: {0} total: {1}  fraction: {2}"
      .format(sgrb2_cloud_mask.sum(),
              nob2cube.mask.include()[0,:,:].sum(),
              sgrb2_cloud_mask.sum() / nob2cube.mask.include()[0,:,:].sum(),
             ))

sp2 = pyspeckit.Spectrum(xarr=full_cubes.cube_merge_high.spectral_axis,
                         data=spd_nob2,
                         xarrkwargs={'velocity_convention': 'radio'},
                        )
sp2.xarr.refX = full_cubes.pcube_merge_high.xarr.refX
sp2.xarr.refX_unit = full_cubes.pcube_merge_high.xarr.refX_unit
sp2.xarr.convert_to_unit(u.GHz)
sp2.specname = 'Whole CMZ'
sp2.unit = 'K'


sp2.plotter()
sp2.specfit.Registry.add_fitter('h2co_simple', simple_fitter2, 6,
                                multisingle='multi')
sp2.baseline(exclude=[151, 761, 916, 1265], selectregion=True,
             highlight_fitregion=True, xtype='pixel', subtract=True)

sp2.error[:] = sp2.data[761:916].std()

sp2.plotter()
sp2.specfit(fittype='h2co_simple',
            guesses=[0.06, 10, 20, 0.5, 0.7, 0.03],
            limited=[(True,True)] * 6,
            limits=[(0,20),[-150,150],(1, 60),(0,1),(0.3,1.1),(0,1e5)],
           )

sp2.specfit.plot_fit(show_components=True)
sp2.plotter.savefig(paths.fpath('simple/WholeCMZ_NoSgrB2_6parameter.pdf'),
                    bbox_inches='tight')

b2_to_total = sgrb2_cloud_mask.sum() / full_cubes.cube_merge_high.with_mask(~sgrb2_cloud_mask).mask[0,:,:].include().sum()

spd_b2 = full_cubes.cube_merge_high.with_mask(sgrb2_cloud_mask).mean(axis=(1,2)).value
sp_b2 = pyspeckit.Spectrum(xarr=full_cubes.cube_merge_high.spectral_axis,
                           data=spd_b2,
                           xarrkwargs={'velocity_convention': 'radio'},
                          )
sp_b2.xarr.refX = full_cubes.pcube_merge_high.xarr.refX
sp_b2.xarr.refX_unit = full_cubes.pcube_merge_high.xarr.refX_unit
sp_b2.xarr.convert_to_unit(u.GHz)
sp_b2.specname = 'Sgr B2'
sp_b2.unit = 'K'


sp_b2_scaled = sp_b2*b2_to_total
ylim = sp2.plotter.axis.get_ylim()
sp_b2_scaled.plotter(axis=sp2.plotter.axis, color='r', clear=False)
sp2.plotter.axis.set_ylim(*ylim)

sp2.plotter.axis.set_title("Whole CMZ and Sgr B2")

sp2.plotter.savefig(paths.fpath('simple/WholeCMZ_noSgrB2_withSgrB2_6parameter.pdf'),
                    bbox_inches='tight')


b2_to_total2 = sgrb2_cloud_mask.sum() / full_cubes.cube_merge_high.mask[0,:,:].include().sum()
sp_b2_scaled2 = sp_b2*b2_to_total2
ylim = sp.plotter.axis.get_ylim()
sp_b2_scaled.plotter(axis=sp.plotter.axis, color='r', clear=False)
sp.plotter.axis.set_ylim(*ylim)
sp.plotter.axis.set_title("Whole CMZ and Sgr B2")
sp.plotter.savefig(paths.fpath('simple/WholeCMZ_withSgrB2_6parameter.pdf'),
                   bbox_inches='tight')

sgrb2_N_mask = ((lon-0.6771*u.deg)**2 + (lat+0.0274*u.deg)**2)**0.5 < 0.5*u.arcmin
spd_b2_N = full_cubes.cube_merge_high.with_mask(sgrb2_N_mask).mean(axis=(1,2)).value
sp_b2_N = pyspeckit.Spectrum(xarr=full_cubes.cube_merge_high.spectral_axis,
                             data=spd_b2_N,
                             xarrkwargs={'velocity_convention': 'radio'},
                            )
sp_b2_N.xarr.refX = full_cubes.pcube_merge_high.xarr.refX
sp_b2_N.xarr.refX_unit = full_cubes.pcube_merge_high.xarr.refX_unit
sp_b2_N.xarr.convert_to_unit(u.GHz)
sp_b2_N.specname = 'Sgr B2 N'
sp_b2_N.unit = 'K'

sgrb2_M_mask = ((lon-0.6668*u.deg)**2 + (lat+0.03508*u.deg)**2)**0.5 < 0.5*u.arcmin
spd_b2_M = full_cubes.cube_merge_high.with_mask(sgrb2_M_mask).mean(axis=(1,2)).value
sp_b2_M = pyspeckit.Spectrum(xarr=full_cubes.cube_merge_high.spectral_axis,
                             data=spd_b2_M,
                             xarrkwargs={'velocity_convention': 'radio'},
                            )
sp_b2_M.xarr.refX = full_cubes.pcube_merge_high.xarr.refX
sp_b2_M.xarr.refX_unit = full_cubes.pcube_merge_high.xarr.refX_unit
sp_b2_M.xarr.convert_to_unit(u.GHz)
sp_b2_M.specname = 'Sgr B2 M'
sp_b2_M.unit = 'K'

N_to_B2 = sgrb2_N_mask.sum() / sgrb2_cloud_mask.sum()
sp_b2_N_scaled = sp_b2_N * N_to_B2

fig4 = pl.figure(4)
fig4.clf()
sp_b2.plotter(axis=fig4.gca(), color='red')
ylim = sp_b2.plotter.axis.get_ylim()
sp_b2_N_scaled.plotter(axis=fig4.gca(), clear=False, color='b')
sp_b2.plotter.axis.set_ylim(*ylim)

sp_b2.plotter.axis.set_title("Sgr B2 and Sgr B2 N")
sp_b2.plotter.savefig(paths.fpath('simple/SgrB2_with_SgrB2N.pdf'),
                      bbox_inches='tight')

fig5 = pl.figure(5)
fig5.clf()
sp_b2_N.plotter(axis=fig5.gca(), color='b')
sp_b2_M.plotter(axis=sp_b2_N.plotter.axis, clear=False, color='g')
sp_b2_N.plotter.axis.set_title("Sgr B2 M and N")
sp_b2_N.plotter.savefig(paths.fpath('simple/SgrB2MandN.pdf'),
                        bbox_inches='tight')


sp_b2.plotter()
sp_b2.specfit.Registry.add_fitter('h2co_simple', simple_fitter2, 6,
                                  multisingle='multi')
sp_b2.baseline(exclude=[151, 761, 916, 1265], selectregion=True,
               highlight_fitregion=True, xtype='pixel', subtract=True)

sp_b2.error[:] = sp_b2.data[761:916].std()

sp_b2.plotter()
sp_b2.specfit(fittype='h2co_simple',
              guesses=[0.06, 10, 20, 0.5, 0.7, 0.03],
              limited=[(True,True)] * 6,
              limits=[(0,20),[-150,150],(1, 60),(0,1),(0.3,1.1),(0,1e5)],
             )

sp_b2.specfit.plot_fit(show_components=True)
sp_b2.plotter.savefig(paths.fpath('simple/SgrB2_6parameter.pdf'),
                      bbox_inches='tight')

from constrain_parameters import paraH2COmodel

mf = paraH2COmodel()

r321303 = sp.specfit.parinfo.values[3]
er321303 = sp.specfit.parinfo.errors[3]
mf.set_constraints(ratio321303=r321303, eratio321303=er321303,
                   logh2column=22,
                   elogh2column=1,
                   logabundance=np.log10(1.2e-9),
                   elogabundance=1, mindens=4, linewidth=10,
                   taline303=sp.specfit.parinfo.values[0],
                   etaline303=sp.specfit.parinfo.errors[0],
                   taline321=sp.specfit.parinfo.values[0]*r321303,
                   etaline321=sp.specfit.parinfo.errors[0],
                   )
constraints = mf.get_parconstraints()
import pylab as pl
pl.gcf().clf()
mf.parplot1d_all()
pl.savefig(paths.fpath('simple/WholeCMZ_1dParConstraints.pdf'),
                   bbox_inches='tight')
