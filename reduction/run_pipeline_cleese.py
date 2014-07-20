import make_apex_cubes
from os.path import join
from astropy import log
import time
import subprocess

root = '/scratch/aginsbur/apex/'
rawpath = join(root,'raw/')
reducedpath = join(root,'reduced/')
make_apex_cubes.june2013datapath = rawpath
make_apex_cubes.june2013path = join(reducedpath,'june2013/')
#make_apex_cubes.april2014path = join(reducedpath,'april2014/')
make_apex_cubes.april2014path = raw
make_apex_cubes.h2copath = join(reducedpath, 'h2co_cubes/')
make_apex_cubes.mergepath = join(reducedpath, 'merged_datasets/')
make_apex_cubes.aorawpath = rawpath
make_apex_cubes.aopath = join(reducedpath, '2010_reduced/')
make_apex_cubes.diagplotdir = join(root,'diagnostic_plots/')

try:
    label = "_v"+subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip()
except CalledProcessError:
    label = ""

logfile = ("".join([time.strftime("apexcmzpipeline{0}_%y_%m_%d_%H:%M:%S"),".log"])).format(label)

with log.log_to_file(logfile):
    make_apex_cubes.do_everything()


#/scratch/aginsbur/apex/reduced/april2014/M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-08.apex
#/scratch/aginsbur/apex/APEX_CMZ_H2CO/raw_data/M-093.F-0009-2014-2014-05-08.apex
#/scratch/aginsbur/apex/raw/M-093.F-0009-2014-2014-05-08.apex
#/scratch/aginsbur/apex/raw/M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-08.apex
