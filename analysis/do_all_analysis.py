execfile("make_apex_cubes.py")
do_postprocessing()

execfile("pyradex_h2comm_grid.py")

execfile("redo_dendro.py")
# execfile("dendro_mask.py")
# make_dend_303()
# execfile("dendro_temperature.py")
# do_dendro_temperatures_both()
# execfile("dendrotem_plots.py")
# execfile("make_piecewise_temcube.py")

execfile("make_ratiotem_cubesims.py")
execfile("make_ratio_integ.py")
make_ratio_integ()
make_ratio_max()
execfile("individual_spectra.py")
execfile("param_plots_densgt4.py")
