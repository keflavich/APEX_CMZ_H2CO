import numpy as np
import time
import warnings
from astropy import log
from astrodendro import Dendrogram,ppv_catalog
from astropy.table import Table, Column
from paths import hpath,tpath
from masked_cubes import cube303,cube303sm,cube321,cube321sm
import flag_other_lines
from astropy.utils.console import ProgressBar

if 'dend' not in locals():
    t0 = time.time()
    log.info("Loading dendrogram from file.")
    dend = Dendrogram.load_from(hpath("DendroMask_H2CO303202.hdf5"))
    log.info("Loaded dendrogram from file in {0:0.1f} seconds.".format(time.time()-t0))
    dend.wcs = cube303.wcs

if 'dendsm' not in locals():
    t0 = time.time()
    log.info("Loading dendrogram from file.")
    dendsm = Dendrogram.load_from(hpath("DendroMask_H2CO303202_smooth.hdf5"))
    log.info("Loaded dendrogram from file in {0:0.1f} seconds.".format(time.time()-t0))
    dendsm.wcs = cube303sm.wcs

# Removed: the 321 signal extraction approach never worked nicely
# if 'dend321' not in locals():
#     t0 = time.time()
#     log.info("Loading dendrogram from file.")
#     dend321 = Dendrogram.load_from(hpath("DendroMask_H2CO321220.hdf5"))
#     log.info("Loaded dendrogram from file in {0:0.1f} seconds.".format(time.time()-t0))
#     dend321.wcs = cube321.wcs
# 
# if 'dend321sm' not in locals():
#     t0 = time.time()
#     log.info("Loading dendrogram from file.")
#     dend321sm = Dendrogram.load_from(hpath("DendroMask_H2CO321220sm.hdf5"))
#     log.info("Loaded dendrogram from file in {0:0.1f} seconds.".format(time.time()-t0))
#     dend321sm.wcs = cube321sm.wcs

# Workaround / bugfix for https://github.com/astropy/astropy/issues/2974
def fix_logical(t):
    newcols = []
    for col in t.columns.values():
        if col.dtype.str.endswith('S5'):
            falses = col == 'False'
            trues = col == 'True'
            if np.all(falses | trues):
                col = t.ColumnClass(trues, name=col.name)
        newcols.append(col)
    return Table(newcols)

dend_sm = dendsm

try:
    catalog = fix_logical(Table.read(tpath('PPV_H2CO_Temperature.ipac'),
                                     format='ascii.ipac', guess=False))
    catalog_sm = fix_logical(Table.read(tpath('PPV_H2CO_Temperature_smooth.ipac'),
                             format='ascii.ipac', guess=False))
    catalogsm = catalog_sm

except IOError:
    warnings.warn("Need to run dendro_temperature:do_dendro_temperatures_both"
                  " to get the PPV H2CO catalogs first")


#try:
#    catalog321 = Table.read(hpath('PPV_H2CO_Temperature_321selected.ipac'),
#                            format='ascii.ipac', guess=False)
#    catalog321_sm = Table.read(hpath('PPV_H2CO_Temperature_321selected_smooth.ipac'),
#                               format='ascii.ipac', guess=False)
#    catalog321sm = catalog321_sm
#
#except IOError:
#    warnings.warn("Need to run dendro_temperature:do_321_dendro_temperatures_both"
#                  " to get the PPV H2CO catalogs first")

# Mark all non-H2CO lines as 'bad'
# This is done manually since there is overlap between HC3N and H2CO 303 in
# velocity
# Also mark Sgr B2's absorption/emission as 'bad'
flag_other_lines.flag_hc3n(dend, catalog=catalog, smooth=False)
flag_other_lines.flag_hc3n(dendsm, catalog=catalog_sm, smooth=True)
flag_other_lines.flag_absorption(dend, catalog=catalog, smooth=False)
flag_other_lines.flag_absorption(dendsm, catalog=catalog_sm, smooth=True)

if 'DespoticTem' not in catalog.colnames:
    from despotic_heating import tkin_all
    print("Adding DESPOTIC-derived temperatures to dendrograms.")
    from astropy import units as u
    import gaussian_correction
    gcorfactor = gaussian_correction.gaussian_correction(catalog['Smin303']/catalog['Smax303'])
    # use 2*reff because we care about diameter, not radius
    # then deconvolve out the beam
    # (see paper: 28" -> 30" beam because of gridding_
    beam_pc = (30/(np.sqrt(8*np.log(2)))*u.arcsec*8.5*u.kpc).to(u.pc, u.dimensionless_angles())
    dtems = [tkin_all(density=10**row['density_chi2']*u.cm**-3,
                      sigma=row['v_rms']*u.km/u.s*gf,
                      lengthscale=2*((row['reff']*u.pc*gf)**2-(beam_pc)**2)**0.5,
                      gradient=5*u.km/u.s/u.pc, #min(5,row['v_rms']/row['reff'])*u.km/u.s/u.pc,
                      tdust=row['higaldusttem']*u.K,
                      crir=1e-17*u.s**-1,
                      ISRF=1,
                      tdust_rad=(row['higaldusttem']*u.K *
                                 (1-np.exp(-(10**row['logh2column']/1e24)))))
             for row,gf in ProgressBar(zip(catalog, gcorfactor))]
    catalog.add_column(Column(name='DespoticTem', data=dtems))
    catalog.add_column(Column(name='gausscorrfactor', data=gcorfactor))

    catalog.write(tpath('PPV_H2CO_Temperature.ipac'), format='ascii.ipac')

