import os

june2013datapath = '/Users/adam/work/h2co/apex/june2013/raw/M-091.F-0019-2013/'
june2013path = '/Users/adam/work/h2co/apex/june2013/'
april2014path = '/Users/adam/work/h2co/apex/april2014/'
h2copath = '/Users/adam/work/h2co/apex/h2co_cubes/'
mergepath = '/Users/adam/work/h2co/apex/merged_datasets/'
aorawpath = '/Users/adam/work/h2co/apex/2010_reduced/2010_raw/'
aopath = '/Users/adam/work/h2co/apex/2010_reduced/'
diagplotdir = '/Users/adam/work/h2co/apex/diagnostic_plots/'
figurepath = '/Users/adam/work/apex_cmz_h2co/tex/figures/'
regpath = '/Users/adam/work/apex_cmz_h2co/regions/'
gridpath = '/Users/adam/work/h2co/radex/thermom/'

def mpath(x, mergepath=mergepath):
    return os.path.join(mergepath,x)

def gpath(x, gridpath=gridpath):
    return os.path.join(gridpath, x)

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)
