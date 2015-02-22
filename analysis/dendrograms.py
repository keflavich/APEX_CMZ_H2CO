import time
import warnings
from astropy import log
from astrodendro import Dendrogram,ppv_catalog
from astropy.table import Table
from paths import hpath
from masked_cubes import cube303,cube303sm,cube321,cube321sm
import flag_other_lines

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


dend_sm = dendsm

try:
    catalog = Table.read(hpath('PPV_H2CO_Temperature.ipac'), format='ascii.ipac',
                         guess=False)
    catalog_sm = Table.read(hpath('PPV_H2CO_Temperature_smooth.ipac'), format='ascii.ipac',
                         guess=False)
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
