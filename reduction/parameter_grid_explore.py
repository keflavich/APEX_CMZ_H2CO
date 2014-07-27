import socket
import itertools
import inspect
import os
import time
import operator
import spectral_cube
import numpy as np
from astropy import units as u
from astropy.io import fits
import pylab as pl

dirpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 

if 'eso-macbook' in socket.gethostname():
    execfile(os.path.join(dirpath,'run_pipeline_cyg.py'))
    outpath = '/Volumes/passport/apex/parameter_tests/'
elif 'cleese' in socket.gethostname():
    execfile(os.path.join(dirpath,'run_pipeline_cleese.py'))
    outpath = '/scratch/aginsbur/apex/parameter_tests/'
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

short_names = {'timewise_pca':'tp',
               'pca_clean':'pca',
               'scanblsub':'blsub',
               'subspectralmeans':'msub',
               'subtract_time_average':'tsub',
               'automask':'auto'}

def build_cube_function_2014(test_dataset, sourcename, outpath, pars, suffix):
    make_apex_cubes.build_cube_2014(source_name,
                                    datasets=[test_dataset],
                                    lowhigh='low',
                                    outpath=outpath,
                                    datapath=make_apex_cubes.april2014path,
                                    extra_suffix="_"+suffix,
                                    **pars)

def parameter_grid_explore_2014():
    return parameter_grid_explore(build_cube_function=build_cube_function_2014,
                                  test_dataset='M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-12',
                                  source_name='MAP_055',
                                  fname_template='APEX_H2CO_2014_{source_name}_{lowhigh}_{suffix}_sub.fits',
                                  lowhigh='low',
                                  outpath=outpath,)


def build_cube_function_2013(test_dataset, sourcename, outpath, pars, suffix):
    make_apex_cubes.build_cube_2013(datasets=[test_dataset],
                                    lowhigh='high',
                                    outpath=outpath,
                                    datapath=make_apex_cubes.june2013datapath,
                                    extra_suffix="_"+suffix,
                                    **pars)

def parameter_grid_explore_2013():
    return parameter_grid_explore(build_cube_function=build_cube_function_2013,
                                  test_dataset='M-091.F-0019-2013-2013-06-08',
                                  source_name=None,
                                  lowhigh='high',
                                  fname_template='APEX_H2CO_2013_{lowhigh}_{suffix}_sub.fits',
                                  outpath=outpath,)


def build_cube_function_ao(test_dataset, sourcename, outpath, pars, suffix):
    make_apex_cubes.build_cube_ao(window='high',
                                  datasets=[test_dataset],
                                  outpath=outpath,
                                  datapath=make_apex_cubes.aorawpath,
                                  extra_suffix="_"+suffix,
                                  **pars)

def parameter_grid_explore_ao():
    return parameter_grid_explore(build_cube_function=build_cube_function_ao,
                                  test_dataset='O-085.F-9311A-2010',
                                  source_name=None,
                                  lowhigh='high',
                                  fname_template='APEX_H2CO_Ao_{lowhigh}_{suffix}_sub.fits',
                                  outpath=outpath,)




def parameter_grid_explore(build_cube_function,
                           test_dataset='M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-12',
                           source_name='MAP_055', 
                           fname_template='',
                           lowhigh='',
                           outpath='/Volumes/passport/apex/parameter_tests/',):

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    results = {}

    for pars in par_sets_pruned:
        suffix = "_".join(short_names[k]+str(int(v)) for k,v in pars.iteritems()) 
        log.info(suffix)
        t0 = time.time()
        build_cube_function(test_dataset=test_dataset,
                            sourcename=source_name,
                            outpath=outpath,
                            pars=pars,
                            suffix=suffix)

        dt = time.time() - t0
        log.info("Finished {0} in {1}s".format(suffix, dt))

        fpath = os.path.join(outpath,
                             fname_template.format(source_name=source_name,
                                                   lowhigh=lowhigh,
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

    times = {r:results[r]['time'] for r in results}
    # http://stackoverflow.com/questions/613183/sort-a-python-dictionary-by-value
    print sorted(times.iteritems(), key=operator.itemgetter(1))

    return results

def plot_pargrid_results(results, outpath=outpath):
    for ii in range(1,7):
        pl.close(ii)
        pl.figure(ii, figsize=(20,20))
        pl.clf()

    for ii,suffix in enumerate(results):

        title = suffix.replace("abcdefghijklmnopqrstuvwxyz_","")
        pl.figure(1)
        pl.subplot(6,5,ii+1)
        val = results[suffix]['integ'].value
        valok = np.isfinite(val)
        pl.imshow(val,
                  vmin=np.percentile(val[valok],1),
                  vmax=np.percentile(val[valok],99),)
        pl.title(title)
        pl.xticks([])
        pl.yticks([])
        pl.suptitle("Integrated")

        pl.figure(2)
        pl.subplot(6,5,ii+1)
        val = results[suffix]['peak'].data
        pl.imshow(val,
                  vmin=np.percentile(val[valok],1),
                  vmax=np.percentile(val[valok],99),)
        pl.title(title)
        pl.xticks([])
        pl.yticks([])
        pl.suptitle("Peak")

        pl.figure(3)
        pl.subplot(6,5,ii+1)
        pl.plot(results[suffix]['spec'].value)
        #pl.ylim(-0.005,0.007)
        pl.xlim(1500,2500)
        pl.title(title)
        pl.xticks([])
        pl.yticks([])

        pl.figure(4)
        pl.subplot(6,5,ii+1)
        val = results[suffix]['integ'].value-results['blsub0_auto0_msub0_pca0_tsub0_tp0']['integ'].value
        pl.imshow(val,
                  vmin=np.percentile(val[valok],5),
                  vmax=np.percentile(val[valok],95),)
        pl.title(title)
        pl.xticks([])
        pl.yticks([])
        pl.suptitle("Integrated: Difference from Fiducial (no cleaning)")

        pl.figure(5)
        pl.subplot(6,5,ii+1)
        val = results[suffix]['integ'].value-results['blsub0_auto0_msub1_pca1_tsub0_tp1']['integ'].value
        pl.imshow(val,
                  vmin=np.percentile(val[valok],1),
                  vmax=np.percentile(val[valok],99),)
        pl.title(title)
        pl.xticks([])
        pl.yticks([])
        pl.suptitle("Integrated: Difference from TP1")
        
        pl.figure(6)
        pl.subplot(6,5,ii+1)
        pl.plot(results[suffix]['spec'].value-results['blsub0_auto0_msub0_pca0_tsub0_tp0']['spec'].value)
        pl.title(title)
        pl.xticks([])
        pl.yticks([])
        #pl.ylim(-0.005,0.007)
        pl.xlim(1500,2500)
        pl.suptitle("Spectrum: Difference from Fiducial (no cleaning)")


    pl.figure(1)
    pl.savefig(os.path.join(outpath, "IntegratedMaps_ParameterExamination.pdf"), bbox_inches='tight')
    pl.figure(2)
    pl.savefig(os.path.join(outpath, "PeakMaps_ParameterExamination.pdf"), bbox_inches='tight')
    pl.figure(3)
    pl.savefig(os.path.join(outpath, "Spectra_ParameterExamination.pdf"), bbox_inches='tight')
    pl.figure(4)
    pl.savefig(os.path.join(outpath, "IntegratedDiffMaps_Fiducial_ParameterExamination.pdf"), bbox_inches='tight')
    pl.figure(5)
    pl.savefig(os.path.join(outpath, "IntegratedDiffMaps_Tp1_ParameterExamination.pdf"), bbox_inches='tight')
    pl.figure(6)
    pl.savefig(os.path.join(outpath, "SpectrumDiff_Fiducial_ParameterExamination.pdf"), bbox_inches='tight')

    """
    [('blsub0_auto0_msub0_pca1_tsub1_tp0', 630.4537680149078),
     ('blsub0_auto0_msub0_pca1_tsub1_tp1', 685.7508239746094),
     ('blsub1_auto1_msub0_pca1_tsub0_tp1', 717.190927028656),
     ('blsub0_auto0_msub0_pca0_tsub0_tp0', 724.5496709346771),
     ('blsub1_auto1_msub0_pca1_tsub0_tp0', 710.7918629646301),
     ('blsub1_auto0_msub0_pca1_tsub0_tp0', 841.5442039966583),
     ('blsub1_auto0_msub0_pca1_tsub0_tp1', 880.426764011383),
     ('blsub0_auto0_msub1_pca0_tsub1_tp0', 608.3833608627319),
     ('blsub0_auto0_msub0_pca1_tsub0_tp1', 632.175253868103),
     ('blsub0_auto0_msub0_pca1_tsub0_tp0', 660.4072258472443),
     ('blsub1_auto2_msub1_pca0_tsub0_tp0', 825.590637922287),
     ('blsub0_auto0_msub1_pca1_tsub0_tp0', 654.482225894928),
     ('blsub0_auto0_msub1_pca1_tsub0_tp1', 636.7380459308624),
     ('blsub1_auto0_msub1_pca0_tsub0_tp0', 645.0383491516113),
     ('blsub1_auto1_msub1_pca1_tsub0_tp0', 703.2395470142365),
     ('blsub1_auto2_msub1_pca1_tsub0_tp0', 786.808100938797),
     ('blsub1_auto2_msub1_pca1_tsub0_tp1', 1204.1531519889832),
     ('blsub0_auto0_msub1_pca0_tsub0_tp0', 761.8484349250793),
     ('blsub1_auto2_msub0_pca0_tsub0_tp0', 907.2306189537048),
     ('blsub1_auto1_msub1_pca1_tsub0_tp1', 688.1312038898468),
     ('blsub0_auto0_msub0_pca0_tsub1_tp0', 652.8386740684509),
     ('blsub1_auto1_msub0_pca0_tsub0_tp0', 692.613874912262),
     ('blsub1_auto2_msub0_pca1_tsub0_tp1', 1603.2192330360413),
     ('blsub0_auto0_msub1_pca1_tsub1_tp1', 645.9985570907593),
     ('blsub0_auto0_msub1_pca1_tsub1_tp0', 691.3041341304779),
     ('blsub1_auto2_msub0_pca1_tsub0_tp0', 843.9561231136322),
     ('blsub1_auto0_msub1_pca1_tsub0_tp1', 1178.8006839752197),
     ('blsub1_auto0_msub1_pca1_tsub0_tp0', 681.291109085083),
     ('blsub1_auto1_msub1_pca0_tsub0_tp0', 1071.6641550064087),
     ('blsub1_auto0_msub0_pca0_tsub0_tp0', 713.5719909667969)]
    """
