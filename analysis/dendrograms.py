import time
from astropy import log
from astrodendro import Dendrogram,ppv_catalog
from astropy.table import Table
from paths import hpath

if 'dend' not in locals():
    t0 = time.time()
    log.debug("Loading dendrogram from file.")
    dend = Dendrogram.load_from(hpath("DendroMask_H2CO303202_signal_to_noise.hdf5"))
    log.debug("Loaded dendrogram from file in {0:0.1f} seconds.".format(time.time()-t0))

if 'dendsm' not in locals():
    t0 = time.time()
    log.debug("Loading dendrogram from file.")
    dendsm = Dendrogram.load_from(hpath("DendroMask_H2CO303202_smooth_signal_to_noise.hdf5"))
    log.debug("Loaded dendrogram from file in {0:0.1f} seconds.".format(time.time()-t0))

catalog = Table.read(hpath('PPV_H2CO_Temperature.ipac'), format='ascii.ipac',
                     guess=False)
#catalog_sm = Table.read(hpath('PPV_H2CO_Temperature_smooth.ipac'), format='ascii.ipac',
#                     guess=False)
