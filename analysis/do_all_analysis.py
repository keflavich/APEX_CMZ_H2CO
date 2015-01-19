from astropy import log
from datetime import datetime
import time
log.setLevel('INFO')
t0 = time.time()

with log.log_to_file("all_analysis_{0}.log".format(datetime.now().isoformat())):

    log.info("Starting make_apex_cubes postprocessing")
    execfile("make_apex_cubes.py")
    do_postprocessing()
    extract_co_subcubes(mergepath=mergepath)

    log.info("Creating pyradex grid.  dt={0}".format(time.time()-t0))
    execfile("pyradex_h2comm_grid.py")

    log.info("Redo dendro.  dt={0}".format(time.time()-t0))
    execfile("redo_dendro.py")
    # execfile("dendro_mask.py")
    # make_dend_303()
    # execfile("dendro_temperature.py")
    # do_dendro_temperatures_both()
    # execfile("dendrotem_plots.py")
    # execfile("make_piecewise_temcube.py")

    log.info("Make ratiotem cubesims.  dt={0}".format(time.time()-t0))
    execfile("make_ratiotem_cubesims.py")
    log.info("Make ratiotem integ.  dt={0}".format(time.time()-t0))
    execfile("make_ratio_integ.py")
    make_ratio_integ()
    make_ratio_max()

    log.info("Individual spectra.  dt={0}".format(time.time()-t0))
    execfile("individual_spectra.py")
    execfile("param_plots_densgt4.py")

    execfile("orbit_pv.py")
    execfile("backstream_pv.py")
    execfile("wholecmz.py")

    execfile("do_all_plots.py")
