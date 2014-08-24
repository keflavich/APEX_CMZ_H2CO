import pyregion
import numpy as np
import pyspeckit
from astropy import table
import spectral_cube
from paths import h2copath, mergepath, figurepath, regpath, analysispath
import os
from pyspeckit_fitting import simplemodel, simple_fitter, simple_fitter2
try:
    from pyspeckit_fitting import h2co_radex_fitter
    radexfit=True
except ImportError:
    radexfit=False
from astropy.io import fits
import photutils
from astropy import coordinates
from astropy import units as u
from astropy.utils.console import ProgressBar


if 'cube' not in locals():
    # Use the individual spectral line cubes because they have been more
    # cleanly baselined
    # (unfortunately, this doesn't appar to work)
    #h2co303 = pyspeckit.Cube(mpath('APEX_H2CO_303_202_bl.fits'))
    #h2co322 = pyspeckit.Cube(mpath('APEX_H2CO_322_221_bl.fits'))
    #h2co321 = pyspeckit.Cube(mpath('APEX_H2CO_321_220_bl.fits'))
    #cube = pyspeckit.CubeStack([h2co303,h2co321,h2co322])
    #cube.xarr.refX = 218222190000.0
    #cube.xarr.refX_units = 'Hz'
    #cube = pyspeckit.Cube(os.path.join(mergepath, 'APEX_H2CO_merge_high_sub.fits'))
    etamb = 0.75 # http://www.apex-telescope.org/telescope/efficiency/
    #cube.cube /= etamb
    noise = fits.getdata(mergepath+'APEX_H2CO_merge_high_sub_noise.fits') / etamb
    noisehdr = fits.getheader(mergepath+'APEX_H2CO_merge_high_sub_noise.fits')
    #errorcube = noise[None,:,:] * np.ones(cube.cube.shape)

    cube = spectral_cube.SpectralCube.read(os.path.join(mergepath, 'APEX_H2CO_merge_high_sub.fits'))
    cube._data /= etamb

    spectra = {}

# Not necessary:
#cube.Registry.add_fitter('h2co_mm_radex', h2co_radex_fitter, 5,
#                         multisingle='multi')
#cube.Registry.add_fitter('h2co_simple', simple_fitter, 4, multisingle='multi')

#spectra = {
#        'brick1': cube.get_apspec((0.2426,0.0081,30), coordsys='galactic', wunit='arcsec'),
#        'brick2': cube.get_apspec((0.2583,0.0181,30), coordsys='galactic', wunit='arcsec'),
#        '20kms': cube.get_apspec((3.59977500e+02,  -6.25e-02, 30), coordsys='galactic', wunit='arcsec'),
#    }
#
#pars = {
#    'brick1': {'ncomp': 1},
#    'brick2': {'ncomp': 2},
#}
parmap_simple = {'ampH2CO':'AMPLITUDE',
          'ampCH3OH':'AMPCH3OH',
          'width':'WIDTH',
          'center':'VELOCITY',
          'h2coratio':'RATIO',}
parmap_simple2 = {'ampH2CO':'AMPLITUDE',
          'ampCH3OH':'AMPCH3OH',
          'width':'WIDTH',
          'center':'VELOCITY',
          'h2coratio321303':'RATIO321303X',
          'h2coratio322321':'RATIO322321X',}
parmap_radex = {
          'temperature':'TEMPERATURE',
          'density':'DENSITY',
          'column':'COLUMN',
          'denswidth':'WIDTH',
          'denscenter':'CENTER',}

def set_row(parinfo, ncomp, row, parmap):

    for ii in range(ncomp):
        for par in parmap:
            row[par+"_"+str(ii)] = parinfo[parmap[par]+str(ii)].value
            row["e"+par+"_"+str(ii)] = parinfo[parmap[par]+str(ii)].error


regs = pyregion.open(regpath+'spectral_apertures.reg') + pyregion.open(regpath+'target_fields_8x8.reg')
with open(regpath+'spectral_ncomp.txt') as f:
    pars = eval(f.read())

width_min,width_max = 1,15

name_column = table.Column(data=[reg.attr[1]['text'] for reg in regs],
                           name='Source_Name')
lon_column = table.Column(data=[reg.coord_list[0] for reg in regs],
                          name='GLON')
lat_column = table.Column(data=[reg.coord_list[1] for reg in regs],
                          name='GLAT')
columns = [table.Column(name="{ee}{name}_{ii}".format(name=name,
                                                      ii=ii,
                                                      ee=ee),
                        dtype='float',
                        length=len(regs))
           for ii in range(max([pars[p]['ncomp'] for p in pars]))
           for name in ['ampH2CO','ampCH3OH','width','center','h2coratio321303',
                        'h2coratio322321',
                        'density','column','temperature','denscenter','denswidth']
           for ee in ['','e']
          ]
out_table = table.Table([name_column, lon_column, lat_column] + columns)

# TODO: replace this with a permanent path
column_image = fits.open('/Users/adam/work/gc/gcmosaic_column_conv36.fits')[0]
dusttem_image = fits.open('/Users/adam/work/gc/gcmosaic_temp_conv36.fits')[0]

surfdens = []
dusttem = []
log.info("Herschel parameter extraction.")
for reg in ProgressBar(regs):
    mask = pyregion.ShapeList([reg]).get_mask(column_image)
    surfdens.append(column_image.data[mask].mean())
    dusttem.append(dusttem_image.data[mask].mean())

surfdens_column = table.Column(data=surfdens, dtype='float',
                               name='higalcolumndens')
dusttem_column = table.Column(data=dusttem, dtype='float', name='higaldusttem')
out_table.add_column(surfdens_column)
out_table.add_column(dusttem_column)

xarr = pyspeckit.units.SpectroscopicAxis(cube.spectral_axis.value,
                                         unit=str(cube.spectral_axis.unit),
                                         refX=cube.wcs.wcs.restfrq,
                                         refX_units='Hz')


for row_number,reg in enumerate(regs):
    name = reg.attr[1]['text']
    if name not in spectra:
        #sp = cube.get_apspec(reg.coord_list,coordsys='galactic',wunit='degree')
        shape = pyregion.ShapeList([reg])
        mask = shape.get_mask(header=noisehdr, shape=noise.shape)
        scube = cube.subcube_from_ds9region(shape)
        data = scube.apply_numpy_function(np.nanmean, axis=(1,2))
        sp = pyspeckit.Spectrum(data=data, xarr=xarr, header=cube.wcs.to_header())
        sp.specname = reg.attr[1]['text']
        # Error is already computed above; this is an old hack
        #sp.error[:] = sp.stats((218e9,218.1e9))['std']
        spectra[name] = sp
    else:
        sp = spectra[name]


    sp.plotter()

    ncomp = pars[sp.specname]['ncomp']
    velos = pars[sp.specname]['velo']
    spname = sp.specname.replace(" ","_")

    sp.specfit.Registry.add_fitter('h2co_simple', simple_fitter2, 6, multisingle='multi')
    guesses_simple = [x for ii in range(ncomp) 
                      for x in (1,velos[ii],5,0.5,1.0,1)]
    sp.specfit(fittype='h2co_simple', multifit=True,
               guesses=guesses_simple)

    set_row(sp.specfit.parinfo, ncomp, out_table[row_number], parmap=parmap_simple2)


    sp.plotter()
    sp.specfit.plot_fit()
    sp.plotter.savefig(os.path.join(figurepath,
                                    "{0}_fit_4_lines_simple.pdf".format(spname)))

    if radexfit:
        guesses = [x for ii in range(ncomp)
                   for x in (100,14,4.5,
                             sp.specfit.parinfo['VELOCITY{0}'.format(ii)].value,
                             (sp.specfit.parinfo['WIDTH{0}'.format(ii)].value
                              if
                              (sp.specfit.parinfo['WIDTH{0}'.format(ii)].value
                               < width_max and
                               sp.specfit.parinfo['WIDTH{0}'.format(ii)].value
                               > width_min)
                              else 5))
                  ]

        sp.specfit.Registry.add_fitter('h2co_mm_radex', h2co_radex_fitter, 5,
                                 multisingle='multi')
        sp.specfit(fittype='h2co_mm_radex', multifit=True,
                   guesses=guesses,
                   limits=[(10,300),(11,15),(3,5.5),(-105,125),(width_min,width_max)]*ncomp,
                   limited=[(True,True)]*5*ncomp,
                   fixed=[False,False,False,True,True]*ncomp,
                   quiet=False,)
        sp.plotter.savefig(os.path.join(figurepath,
                                        "{0}_fit_h2co_mm_radex.pdf".format(spname)))

        set_row(sp.specfit.parinfo, ncomp, out_table[row_number], parmap=parmap_radex)

    individual_fits=False
    if individual_fits:
        flux = 3.6 # Jy
        col_per_jy = 2e22 # cm^-2
        dvdpc = 5.0 # km/s/pc
        logX = -8.3
        logcol = np.log10(flux*col_per_jy/dvdpc) + logX

        spectra['WarmSpot'].specfit(fittype='h2co_mm_radex', multifit=True,
                                    guesses=[100,logcol,4.5,35,3.0],
                                    limits=[(20,200),(11,15),(3,5.5),(-105,105),(1,5)],
                                    limited=[(True,True)]*5,
                                    fixed=[False,True,True,False,False],
                                    quiet=False,)

        spectra['WarmSpot'].specfit(fittype='h2co_mm_radex', multifit=True,
                                    guesses=[100,logcol+np.log10(2/3.),4.5,27,3.0]+[100,logcol+np.log10(1.0/3.),4.5,53,3.0],
                                    limits=[(20,200),(11,15),(3,5.5),(-105,105),(1,8)]+[(20,200),(11,15),(3,5.5),(-105,105),(1,6)],
                                    limited=[(True,True)]*10,
                                    fixed=[False,True,True,False,False]*2,
                                    quiet=False,)

        flux = 5.0
        logX = -8.5
        logcol = np.log10(flux*col_per_jy/dvdpc / 2.) + logX

        spectra['Brick SW'].specfit(fittype='h2co_mm_radex', multifit=True,
                                    guesses=[133,12.94,5.977,37.17,9.88],
                                    limits=[(20,200),(11,15),(3,6.5),(-105,105),(1,15)],
                                    limited=[(True,True)]*5,
                                    #fixed=[False,False,True,False,False],
                                    fixed=[False,True,True,False,False])

        spectra['50kmsColdExtension'].specfit(fittype='h2co_mm_radex', multifit=True,
                                    guesses=[33,12.94,4.0,24.4,7],
                                    limits=[(20,200),(11,15),(3,6.5),(-105,105),(1,15)],
                                    limited=[(True,True)]*5,
                                    #fixed=[False,False,True,False,False],
                                    fixed=[False,False,True,False,False])

        #spectra['Sgr B2 SW'].specfit(fittype='h2co_mm_radex', multifit=True,
        #                            guesses=[125,14.14,4.0,47.72,5.66]+[125,14.14,4.0,55.72,5.66],
        #                            limits=[(20,200),(11,15),(3,6.5),(-105,105),(1,15)]*2,
        #                            limited=[(True,True)]*5*2,
        #                            #fixed=[False,False,True,False,False],
        #                            fixed=[False,False,True,True,True]+[False,False,True,False,False])



    dopymc = False
    if dopymc:
        import agpy

        sp = spectra['20 kms']
        # SHOULD BE 
        spmc = sp.specfit.get_pymc(use_fitted_values=True, db='hdf5', dbname='h2co_mm_fit_20kmsCld.hdf5')
        #spmc = sp.specfit.fitter.get_pymc(sp.xarr, sp.data, sp.error,
        #                                  use_fitted_values=True, db='hdf5',
        #                                  dbname='h2co_mm_fit_20kmsCld.hdf5')
        spmc.sample(100000)
        agpy.pymc_plotting.hist2d(spmc, 'TEMPERATURE0', 'DENSITY0', doerrellipse=False, clear=True, bins=50, fignum=4)
        agpy.pymc_plotting.hist2d(spmc, 'TEMPERATURE0', 'COLUMN0', doerrellipse=False, clear=True, bins=50,fignum=5)
        agpy.pymc_plotting.hist2d(spmc, 'DENSITY0', 'COLUMN0', doerrellipse=False, clear=True, bins=50,fignum=6)


        pars = dict([(k,spmc.trace(k)[-50:]) for k in sp.specfit.parinfo.keys()])
        sp.plotter.autorefresh=False
        for ii in xrange(0,50):
            sp.specfit.plot_model([pars[k][ii] for k in sp.specfit.parinfo.keys()],
                                  clear=False,
                                  composite_fit_color='r',
                                  plotkwargs={'alpha':0.01})

        sp.plotter.refresh()

out_table.write(os.path.join(analysispath,"fitted_line_parameters.ipac"), format='ascii.ipac')
