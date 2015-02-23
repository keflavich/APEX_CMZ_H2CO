import os
import socket

if socket.gethostname() == 'cleese':
    root = '/scratch/aginsbur/apex/'
else:
    root = os.path.expanduser('~/work')

june2013datapath = os.path.join(root, 'h2co/apex/june2013/raw/M-091.F-0019-2013/')
june2013path = os.path.join(root, 'h2co/apex/june2013/')
april2014path = os.path.join(root, 'h2co/apex/april2014/')
h2copath = os.path.join(root, 'h2co/apex/h2co_cubes/')
mergepath = os.path.join(root, 'h2co/apex/merged_datasets/')
molcubepath = os.path.join(mergepath, 'molecule_cubes/')
aorawpath = os.path.join(root, 'h2co/apex/2010_reduced/2010_raw/')
aopath = os.path.join(root, 'h2co/apex/2010_reduced/')
diagplotdir = os.path.join(root, 'h2co/apex/diagnostic_plots/')
figurepath = os.path.join(root, 'apex_cmz_h2co/tex/figures/')
regpath = os.path.join(root, 'apex_cmz_h2co/regions/')
gridpath = os.path.join(root, 'h2co/radex/thermom/')
analysispath = os.path.join(root, 'apex_cmz_h2co/analysis/')
plotcodepath = os.path.join(root, 'apex_cmz_h2co/plot_codes/')
observingpath = os.path.join(root, 'apex_cmz_h2co/observing/')
tablepath = os.path.join(root, 'apex_cmz_h2co/tables/')

def mpath(x, mergepath=mergepath):
    return os.path.join(mergepath,x)

def gpath(x, gridpath=gridpath):
    return os.path.join(gridpath, x)

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)

def rpath(x, regpath=regpath):
    return os.path.join(regpath, x)

def opath(x, observingpath=observingpath):
    return os.path.join(observingpath, x)

def hpath(x, h2copath=h2copath):
    return os.path.join(h2copath, x)
datapath = h2copath
dpath = hpath

def pcpath(x, plotcodepath=plotcodepath):
    return os.path.join(plotcodepath, x)

def apath(x, analysispath=analysispath):
    return os.path.join(analysispath, x)

def molpath(x, molcubepath=molcubepath):
    return os.path.join(molcubepath,x)

def tpath(x, tablepath=tablepath):
    return os.path.join(tablepath, x)
