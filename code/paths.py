import os

root = os.path.expanduser('~/work')

june2013datapath = os.path.join(root, 'h2co/apex/june2013/raw/M-091.F-0019-2013/')
june2013path = os.path.join(root, 'h2co/apex/june2013/')
april2014path = os.path.join(root, 'h2co/apex/april2014/')
h2copath = os.path.join(root, 'h2co/apex/h2co_cubes/')
mergepath = os.path.join(root, 'h2co/apex/merged_datasets/')
aorawpath = os.path.join(root, 'h2co/apex/2010_reduced/2010_raw/')
aopath = os.path.join(root, 'h2co/apex/2010_reduced/')
diagplotdir = os.path.join(root, 'h2co/apex/diagnostic_plots/')
figurepath = os.path.join(root, 'apex_cmz_h2co/tex/figures/')
regpath = os.path.join(root, 'apex_cmz_h2co/regions/')
gridpath = os.path.join(root, 'h2co/radex/thermom/')
analysispath = os.path.join(root, 'apex_cmz_h2co/analysis/')

def mpath(x, mergepath=mergepath):
    return os.path.join(mergepath,x)

def gpath(x, gridpath=gridpath):
    return os.path.join(gridpath, x)

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)
