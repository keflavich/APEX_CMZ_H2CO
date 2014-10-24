from shfi_otf_pipeline import make_apex_cubes
from os.path import join
from astropy import log
import time
import subprocess
import paths

root = '/scratch/aginsbur/apex/'
paths.root = root
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
make_apex_cubes.molpath = os.path.join(make_apex_cubes.mergepath,
                                       'molecule_cubes/')

try:
    label = "_v"+subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip()
except CalledProcessError:
    label = ""

logfile = ("".join([time.strftime("apexcmzpipeline{0}_%y_%m_%d_%H:%M:%S"),".log"])).format(label)

def do_everything():
    make_apex_cubes.do_everything(mergepath=make_apex_cubes.mergepath,
                                  h2copath=make_apex_cubes.h2copath,
                                  molpath=make_apex_cubes.molpath)
def do_postprocessing():
    make_apex_cubes.do_postprocessing(mergepath=make_apex_cubes.mergepath,
                                      h2copath=make_apex_cubes.h2copath,
                                      molpath=make_apex_cubes.molpath)

#with log.log_to_file(logfile.replace(".log","_do_everything_default.log")):
#   make_apex_cubes.do_everything()

#with log.log_to_file(logfile.replace(".log","_do_everything_default.log")):
#   make_apex_cubes.do_everything(mergefile2='APEX_H2CO_merge_high_default_August2014')

# with log.log_to_file(logfile.replace(".log","_default.log")):
#    make_apex_cubes.make_high_mergecube(mergefile2='APEX_H2CO_merge_high_default')
# with log.log_to_file(logfile.replace(".log","_nopca.log")):
#    make_apex_cubes.make_high_mergecube(pca_clean=False, mergefile2='APEX_H2CO_merge_high_nopca')
# with log.log_to_file(logfile.replace(".log","_timepca.log")):
#    make_apex_cubes.make_high_mergecube(timewise_pca=True, mergefile2='APEX_H2CO_merge_high_timepca')
# with log.log_to_file(logfile.replace(".log","_unclean.log")):
#    make_apex_cubes.make_high_mergecube(pca_clean={'2014':False, '2013':False, 'ao':False},
#                                        scanblsub={'2014':False, '2013':False, 'ao':False},
#                                        timewise_pca={'2014':False, '2013':False, 'ao':False}, 
#                                        mergefile2='APEX_H2CO_merge_high_nopca_noscanblsub')
# with log.log_to_file(logfile.replace(".log","_only_clean_2014.log")):
#    make_apex_cubes.make_high_mergecube(pca_clean={'2014':True, '2013':False, 'ao':False},
#                                        scanblsub={'2014':False, '2013':False, 'ao':False},
#                                        timewise_pca={'2014':True, '2013':False, 'ao':False}, 
#                                        mergefile2='APEX_H2CO_merge_high_pca2014only')
#

#/scratch/aginsbur/apex/reduced/april2014/M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-08.apex
#/scratch/aginsbur/apex/APEX_CMZ_H2CO/raw_data/M-093.F-0009-2014-2014-05-08.apex
#/scratch/aginsbur/apex/raw/M-093.F-0009-2014-2014-05-08.apex
#/scratch/aginsbur/apex/raw/M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-08.apex
