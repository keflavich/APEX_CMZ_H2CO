from shfi_otf_pipeline import make_apex_cubes
from os.path import join
from astropy import log
import time
import subprocess

import os

if os.path.exists('/Volumes/passport/apex/'):
    log.info("Treating /Volumes/passport/ as the root directory.")

    root = '/Volumes/passport/apex/'
    rawpath = join(root,'raw/')
    reducedpath = join(root,'reduced/')
    make_apex_cubes.june2013datapath = rawpath
    make_apex_cubes.june2013path = join(reducedpath,'june2013/')
    #make_apex_cubes.april2014path = join(reducedpath,'april2014/')
    make_apex_cubes.april2014path = rawpath
    make_apex_cubes.h2copath = join(reducedpath, 'h2co_cubes/')
    make_apex_cubes.mergepath = join(reducedpath, 'merged_datasets/')
    make_apex_cubes.aorawpath = rawpath
    make_apex_cubes.aopath = join(reducedpath, '2010_reduced/')
    make_apex_cubes.diagplotdir = join(root,'diagnostic_plots/')
