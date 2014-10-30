import time
from astropy import log
from astrodendro import Dendrogram,ppv_catalog
from astropy.table import Table
from paths import hpath
from masked_cubes import cube303,cube303sm
import flag_other_lines

if 'dend' not in locals():
    t0 = time.time()
    log.info("Loading dendrogram from file.")
    dend = Dendrogram.load_from(hpath("DendroMask_H2CO303202_signal_to_noise.hdf5"))
    log.info("Loaded dendrogram from file in {0:0.1f} seconds.".format(time.time()-t0))
    dend.wcs = cube303.wcs

if 'dendsm' not in locals():
    t0 = time.time()
    log.info("Loading dendrogram from file.")
    dendsm = Dendrogram.load_from(hpath("DendroMask_H2CO303202_smooth_signal_to_noise.hdf5"))
    log.info("Loaded dendrogram from file in {0:0.1f} seconds.".format(time.time()-t0))
    dendsm.wcs = cube303sm.wcs

dend_sm = dendsm

catalog = Table.read(hpath('PPV_H2CO_Temperature.ipac'), format='ascii.ipac',
                     guess=False)
catalog_sm = Table.read(hpath('PPV_H2CO_Temperature_smooth.ipac'), format='ascii.ipac',
                     guess=False)
catalogsm = catalog_sm

# Mark all non-H2CO lines as 'bad'
# This is done manually since there is overlap between HC3N and H2CO 303 in
# velocity
flag_other_lines.flag_dendro(dend, catalog=catalog, smooth=False)
flag_other_lines.flag_dendro(dendsm, catalog=catalog_sm, smooth=True)
