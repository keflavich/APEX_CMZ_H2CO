import socket
import itertools
import inspect
import os
import time

dirpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 

if 'eso-macbook' in socket.gethostname():
    execfile(os.path.join(dirpath,'run_pipeline_cyg.py'))
elif 'cleese' in socket.gethostname():
    execfile(os.path.join(dirpath,'run_pipeline_cleese.py'))
else:
    raise ValueError("Machine {0} not recognized.".format(socket.gethostname()))

parameters = {'timewise_pca': [True, False],
              'pca_clean': [True, False],
              'scanblsub': [True, False],
              'subspectralmeans': [True, False],
              'subtract_time_average': [True, False],
              'automask': [0, 1, 2],}

par_pairs = [[(k,v) for v in parameters[k]] for k in parameters]
par_sets = [dict(x) for x in itertools.product(*par_pairs)]

def prune(pars):
    if pars['timewise_pca'] and not pars['pca_clean']:
        return False
    if pars['automask'] and not pars['scanblsub']:
        return False
    if pars['scanblsub'] and pars['subtract_time_average']:
        return False

    return True

par_sets_pruned = [p for p in par_sets if prune(p)]

test_dataset='M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-12'
dataset='M-093.F-0009-2014-2014-05-12'
source_name ='MAP_055'

short_names = {'timewise_pca':'tp',
               'pca_clean':'pca',
               'scanblsub':'blsub',
               'subspectralmeans':'msub',
               'subtract_time_average':'tsub',
               'automask':'auto'}

dpath = make_apex_cubes.april2014path
outpath = '/Volumes/passport/apex/parameter_tests/'
if not os.path.exists(outpath):
    os.mkdir(outpath)

results = {}

for pars in par_sets_pruned:
    suffix = "_".join(short_names[k]+str(int(v)) for k,v in pars.iteritems()) 
    log.info(suffix)
    t0 = time.time()
    make_apex_cubes.build_cube_2014(source_name,
                                    datasets=[test_dataset],
                                    lowhigh='low',
                                    outpath=outpath,
                                    extra_suffix="_"+suffix,
                                    **pars)

    dt = time.time() - t0
    log.info("Finished {0} in {1}s".format(suffix, dt))

    fpath = os.path.join(outpath,
                         'APEX_H2CO_2014_{source_name}_{lowhigh}_{suffix}_sub.fits'.
                         format(source_name=source_name, lowhigh='low',
                                suffix=suffix))
    cube = spectral_cube.SpectralCube.read(fpath)

    h2cocube = cube.with_spectral_unit(u.km/u.s, rest_value=218.22219*u.GHz, velocity_convention='radio')
    h2cocubecut = h2cocube.spectral_slab(15*u.km/u.s, 65*u.km/u.s)
    inth2co = h2cocubecut.moment0()
    peakh2co = fits.PrimaryHDU(h2cocubecut.max(axis=0).value, header=inth2co.hdu.header)
    spec = h2cocube[:,20:40,20:40].moment0(axis=2).mean(axis=1)
    
    results[suffix] = {'cube': h2cocube,
                       'cubecut': h2cocubecut,
                       'integ': inth2co,
                       'peak': peakh2co,
                       'spec': spec,
                       'time': dt}

import pylab as pl

for ii,suffix in enumerate(results):

    pl.figure(1)
    pl.subplot(6,5,ii+1)
    pl.imshow(results['suffix']['inth2co'].value)

    pl.figure(2)
    pl.subplot(6,5,ii+1)
    pl.imshow(results['suffix']['peakh2co'].value)

    pl.figure(3)
    pl.subplot(6,5,ii+1)
    pl.plot(results['suffix']['spec'].value)
    pl.ylim(-0.005,0.007)
    pl.xlim(1500,2500)
