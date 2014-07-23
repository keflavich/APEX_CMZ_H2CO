import socket
import inspect, os

dirpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 

if 'eso-macbook' in socket.gethostname():
    execfile(os.path.join(dirpath,'run_pipeline_cyg.py'))
elif 'cleese' in socket.gethostname():
    execfile(os.path.join(dirpath,'run_pipeline_cleese.py'))
else:
    raise ValueError("Machine {0} not recognized.".format(socket.gethostname()))


try:
    label = "_v"+subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip()
except CalledProcessError:
    label = ""

logfile = ("".join([time.strftime("apexcmzpipeline{0}_%y_%m_%d_%H:%M:%S"),".log"])).format(label)

with log.log_to_file(logfile):
    for dataset in make_apex_cubes.datasets_2014:
        mapnames = make_apex_cubes.datasets_2014[dataset]

        for mapname in mapnames:
            make_apex_cubes.build_cube_2014(mapname,
                                            datapath=make_apex_cubes.april2014path,
                                            outpath=make_apex_cubes.april2014path,
                                            lowhigh='low', pca_clean=True,
                                            pcakwargs={}, datasets=[dataset])

    for dataset in make_apex_cubes.datasets_2013:
        make_apex_cubes.build_cube_2013(datapath=make_apex_cubes.june2013path,
                                        outpath=make_apex_cubes.june2013path,
                                        lowhigh='low', pca_clean=True,
                                        pcakwargs={}, datasets=[dataset],
                                        extra_suffix=dataset[-5:])

    for dataset in make_apex_cubes.datasets_ao:
        make_apex_cubes.build_cube_ao(window='high', datasets=[dataset],
                                      datapath=make_apex_cubes.aorawpath,
                                      outpath=make_apex_cubes.aopath,
                                      timewise_pca=True, pca_clean=True,
                                      freq=True,)

