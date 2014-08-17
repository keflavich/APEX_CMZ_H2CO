import os
root = '/Users/adam/work/gc/apex/'
datapath = os.path.join(root, 'h2co_cubes')
mergepath = os.path.join(root, 'merged_datasets') 

def dpath(x, datapath=datapath):
    return os.path.join(datapath, x)

def mpath(x, datapath=mergepath):
    return os.path.join(datapath, x)

