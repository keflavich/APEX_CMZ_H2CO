import pyregion
import numpy as np
import pyspeckit
from astropy import table
import spectral_cube
from paths import h2copath, mergepath, figurepath, regpath, analysispath, mpath, hpath
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
from astropy import log
from astropy.utils.console import ProgressBar
import pylab as pl
from higal_gridded import dusttem_image, column_image
import copy

from full_cubes import cube_merge_high as cube
from noise import noise, noisehdr
noiseokmask = np.isfinite(noise)

etamb = 0.75 # http://www.apex-telescope.org/telescope/efficiency/
cube._data /= etamb
noise /= etamb

with open(regpath+'spectral_ncomp.txt') as f:
    pars = eval(f.read())

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
parmap_simple2_spline = {'spline_ampH2CO':'AMPLITUDE',
                         'spline_ampCH3OH':'AMPCH3OH',
                         'spline_width':'WIDTH',
                         'spline_center':'VELOCITY',
                         'spline_h2coratio321303':'RATIO321303X',
                         'spline_h2coratio322321':'RATIO322321X',}
#parmap_radex = {
#          'temperature':'TEMPERATURE',
#          'density':'DENSITY',
#          'column':'COLUMN',
#          'denswidth':'WIDTH',
#          'denscenter':'CENTER',}

def set_row(parinfo, ncomp, rows, parmap):

    assert ncomp == len(rows)

    for ii,row in enumerate(rows):
        row['ComponentID'] = ii
        for par in parmap:
            row[par] = parinfo[parmap[par]+str(ii)].value
            row["e"+par] = parinfo[parmap[par]+str(ii)].error


font_sizes = {1: 20,
              2: 15,
              3: 11,
              4: 8}

def fit_a_spectrum(sp, radexfit=False, write=True, vlimits=(-105,125)):
    sp.plotter.autorefresh=False
    sp.plotter(figure=1)
    ncomp = pars[sp.specname]['ncomp']
    if ncomp == 0:
        log.info("Skipping {0} - no velocity components detected.".format(ncomp))
        return
    returns = [ncomp]
    velos = pars[sp.specname]['velo']
    spname = sp.specname.replace(" ","_")

    if 'Map' in sp.specname or 'box' in sp.specname:
        width_min,width_max = 1,40
    else:
        width_min,width_max = 1,15


    sp.specfit.Registry.add_fitter('h2co_simple', simple_fitter2, 6,
                                   multisingle='multi')
    guesses_simple = [x for ii in range(ncomp) 
                      for x in (sp.data.max(),velos[ii],5,0.5,1.0,sp.data.max())]

    if not(min(velos) > vlimits[0] and max(velos) < vlimits[1]):
        log.warn("A velocity guess {0} is outside limits {1}." .format(velos,
                                                                       vlimits))
        vlimits = (min(velos)-25, max(velos)+25)
        log.warn("Changing limits to {0}".format(vlimits))

    sp.specfit(fittype='h2co_simple', multifit=True,
               guesses=guesses_simple,
               limited=[(True,True)] * 6,
               limits=[(0,20),vlimits,(width_min,width_max),(0,1),(0.3,1.1),(0,1e5)],
              )
    sp.baseline(excludefit=True, subtract=True, highlight_fitregion=True, order=1)

    sp.plotter(clear=True)
    sp.specfit(fittype='h2co_simple', multifit=True,
               guesses=guesses_simple,
               limited=[(True,True)] * 6,
               limits=[(0,20),vlimits,(width_min,width_max),(0,1),(0.3,1.1),(0,1e5)],
              )

    returns.append(copy.copy(sp.specfit.parinfo))

    err = sp.error.mean()

    sp.plotter()
    sp.specfit.plot_fit(show_components=True)
    sp.specfit.annotate(fontsize=font_sizes[ncomp])
    sp.specfit.plotresiduals(axis=sp.plotter.axis, yoffset=-err*5, clear=False,
                             color='#444444', label=False)
    sp.plotter.axis.set_ylim(sp.plotter.ymin-err*5, sp.plotter.ymax)
    sp.plotter.savefig(os.path.join(figurepath,
                                    "simple/{0}_fit_4_lines_simple.pdf".format(spname)))
    if write:
        sp.write(mpath("spectra/{0}_spectrum.fits".format(spname)))

    # This will mess things up for the radexfit (maybe in a good way) but *cannot*
    # be done after the radexfit
    # Set the spectrum to be the fit residualsa.  The linear baseline has
    # already been subtracted from both the data and the residuals
    linear_baseline = sp.baseline.basespec
    sp.baseline.unsubtract()
    sp.baseline.spectofit = sp.specfit.residuals
    sp.baseline.includemask[:] = True # Select ALL residuals
    sp.baseline.fit(spline=True, order=3, spline_sampling=50)
    spline_baseline = sp.baseline.basespec
    sp.data -= spline_baseline + linear_baseline
    sp.baseline.subtracted = True
    sp.error[:] = sp.stats((218.5e9,218.65e9))['std']
    sp.specfit(fittype='h2co_simple', multifit=True,
               guesses=guesses_simple,
               limited=[(True,True)] * 6,
               limits=[(0,1e5),vlimits,(width_min,width_max),(0,1),(0.3,1.1),(0,1e5)],
              )
    sp.plotter()
    sp.plotter.axis.plot(sp.xarr, spline_baseline+linear_baseline-err*5, color='orange',
                         alpha=0.5, zorder=-1, linewidth=2)
    sp.specfit.plot_fit(show_components=True)
    sp.specfit.annotate(fontsize=font_sizes[ncomp])
    sp.specfit.plotresiduals(axis=sp.plotter.axis, yoffset=-err*5, clear=False,
                             color='#444444', label=False)
    sp.plotter.axis.set_ylim(sp.plotter.ymin-err*5, sp.plotter.ymax)
    sp.plotter.savefig(os.path.join(figurepath,
                                    "simple/{0}_fit_4_lines_simple_splinebaselined.pdf".format(spname)))

    returns.append(copy.copy(sp.specfit.parinfo))

    if write:
        sp.write(mpath("spectra/{0}_spectrum_basesplined.fits".format(spname)))

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
                                        "radex/{0}_fit_h2co_mm_radex.pdf".format(spname)))

        returns.append(copy.copy(sp.specfit.parinfo))

    return returns

# Use the individual spectral line cubes because they have been more
# cleanly baselined
# (unfortunately, this doesn't appar to work)
#h2co303 = pyspeckit.Cube(hpath('APEX_H2CO_303_202_bl.fits'))
#h2co322 = pyspeckit.Cube(hpath('APEX_H2CO_322_221_bl.fits'))
#h2co321 = pyspeckit.Cube(hpath('APEX_H2CO_321_220_bl.fits'))
#cube = pyspeckit.CubeStack([h2co303,h2co321,h2co322])
#cube.xarr.refX = 218222190000.0
#cube.xarr.refX_units = 'Hz'
#cube = pyspeckit.Cube(os.path.join(mergepath, 'APEX_H2CO_merge_high_sub.fits'))
#cube.cube /= etamb
#errorcube = noise[None,:,:] * np.ones(cube.cube.shape)

def load_spectra(regs, cube):

    spectra = {}

    xarr = pyspeckit.units.SpectroscopicAxis(cube.spectral_axis.value,
                                             unit=str(cube.spectral_axis.unit),
                                             refX=cube.wcs.wcs.restfrq,
                                             refX_units='Hz')


    for region_number,reg in enumerate(regs):
        name = reg.attr[1]['text']
        log.info("Loading {0}".format(name))
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
            sp.header['ERROR'] = error
            sp.error[:] = sp.stats((218.5e9,218.65e9))['std']
            sp.specname = reg.attr[1]['text']
            # Error is already computed above; this is an old hack
            #sp.error[:] = sp.stats((218e9,218.1e9))['std']
            spectra[name] = sp
            sp.unit = "$T_{MB}$ [K]"
        else:
            sp = spectra[name]

    return spectra

if __name__ == "__main__":

    pl.ioff()
    pl.close(1)
    pl.figure(1).clf()

    radexfit=False # not super useful...

    regs = (pyregion.open(regpath+'spectral_apertures.reg') +
            pyregion.open(regpath+'target_fields_8x8_gal.reg'))

    name_column = table.Column(data=[reg.attr[1]['text'] 
                                     for reg in regs 
                                     for ii in range(pars[reg.attr[1]['text']]['ncomp'])],
                               name='Source_Name')
    comp_id_column = table.Column(data=[0]*name_column.size, name='ComponentID')
    lon_column = table.Column(data=[reg.coord_list[0]
                                    for reg in regs
                                    for ii in range(pars[reg.attr[1]['text']]['ncomp'])
                                   ],
                              name='GLON')
    lat_column = table.Column(data=[reg.coord_list[1] for reg in regs
                                    for ii in range(pars[reg.attr[1]['text']]['ncomp'])
                                   ],
                              name='GLAT')
    columns = [table.Column(name="{ee}{name}".format(name=name, ee=ee),
                            dtype='float',
                            length=name_column.size)
               for name in (parmap_simple2.keys() +# parmap_radex.keys() +
                            parmap_simple2_spline.keys())
               for ee in ['','e']
              ]
    columns += [table.Column(name="{name}".format(name=name, ee=ee),
                            dtype='float',
                            length=name_column.size)
                for name in ['boxwidth', 'boxheight', 'radius', 'area',
                             'posang'] ]
    out_table = table.Table([name_column, comp_id_column, lon_column, lat_column] +
                            columns)

    surfdens = []
    dusttem = []
    log.info("Herschel parameter extraction.")
    herschelok = np.isfinite(column_image.data) & np.isfinite(dusttem_image.data)
    for reg in ProgressBar(regs):
        mask = pyregion.ShapeList([reg]).get_mask(column_image) & herschelok
        for ii in range(pars[reg.attr[1]['text']]['ncomp']):
            surfdens.append(column_image.data[mask].mean()*1e22)
            dusttem.append(dusttem_image.data[mask].mean())

    surfdens_column = table.Column(data=surfdens, dtype='float',
                                   name='higalcolumndens')
    dusttem_column = table.Column(data=dusttem, dtype='float', name='higaldusttem')
    out_table.add_column(surfdens_column)
    out_table.add_column(dusttem_column)

    row_number = 0

    spectra = load_spectra(regs, cube)

    for region_number,reg in enumerate(regs):
        name = reg.attr[1]['text']
        sp = spectra[name]
        log.info("Fitting {0}".format(name))

        returns = fit_a_spectrum(sp, radexfit=radexfit)

        if returns is None:
            continue
        elif radexfit:
            ncomp, pinf1, pinf2, pinf3 = returns
        else:
            ncomp, pinf1, pinf2 = returns

        if reg.name == 'box':
            out_table[row_number:row_number+ncomp]['boxwidth'] = reg.coord_list[2]
            out_table[row_number:row_number+ncomp]['boxheight'] = reg.coord_list[3]
            out_table[row_number:row_number+ncomp]['area'] = reg.coord_list[2] * reg.coord_list[3]
            out_table[row_number:row_number+ncomp]['posang'] = reg.coord_list[4]
        elif reg.name == 'circle':
            out_table[row_number:row_number+ncomp]['area'] = reg.coord_list[2]**2 * np.pi
            out_table[row_number:row_number+ncomp]['radius'] = reg.coord_list[2]
        else:
            raise ValueError("Unsupported region.  Implement it if you want it.")


        set_row(pinf1, ncomp, out_table[row_number:row_number+ncomp],
                parmap=parmap_simple2)
        set_row(pinf2, ncomp, out_table[row_number:row_number+ncomp],
                parmap=parmap_simple2_spline)
        if radexfit and False: # hard-coded out the radex par mapping above
            set_row(pinf3, ncomp,
                    out_table[row_number:row_number+ncomp], parmap=parmap_radex)

        row_number = row_number + ncomp



    out_table.write(os.path.join(analysispath,"fitted_line_parameters.ipac"),
                    format='ascii.ipac')
