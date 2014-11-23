from paths import rpath,mpath,opath
from make_apex_cubes import all_apexfiles,get_source_tel_line,_is_sci, hdr_to_freq
from pyspeckit.spectrum.readers import read_class
from astropy.table import Table
from astropy import log
from astropy.utils.console import ProgressBar
import numpy as np
import os
import pylab as pl

def tsys_data(plot=False):

    if plot:
        fig1 = pl.figure(1)
        fig2 = pl.figure(2)
        fig1.clf()
        fig2.clf()
        ax1 = fig1.gca()
        ax2 = fig2.gca()

    datadict = {}
    tbldict = {}

    for apex_filename in all_apexfiles:
        log.info(apex_filename)
        cl = read_class.ClassObject(apex_filename)

        sourcereg,line,telescopes = get_source_tel_line(apex_filename)
        sci_sources = [source for source in cl.sources
                       if _is_sci(source, sourcereg)]

        datadict[apex_filename] = {t:[] for t in telescopes}

        for telescope in telescopes:
            log.info('{0}: {1}'.format(apex_filename, telescope))
            selection = [x
                         for source in sci_sources
                         for x in cl.select_spectra(telescope=telescope,
                                                    line=line,
                                                    source=source)]

            spdheader = cl.read_observations(selection, progressbar=True)

            datadict[apex_filename][telescope] = zip(*[(sp.std(), h['TSYS'])
                                                       for sp,h in ProgressBar(spdheader)])

        tbl = Table([datadict[apex_filename][t][ii]
                     for t in telescopes
                     for ii in (0,1)],
                    names=[t+"_"+s 
                           for t in telescopes
                           for s in ('STDDEV','TSYS',)],
                    dtype=['float'
                           for t in telescopes
                           for s in ('STDDEV','TSYS',)
                           ])
        log.info(os.path.basename(apex_filename)+"_tsys.fits")
        tbl.write(os.path.basename(apex_filename)+"_tsys.fits", overwrite=True)
        tbldict[apex_filename] = tbl

        if plot:
            ax1.plot(tbl['{0}_TSYS'.format(telescopes[0])],
                     tbl['{0}_STDDEV'.format(telescopes[0])],
                     ',', alpha=0.8)
            ax1.set_xlabel("TSYS")
            ax1.set_ylabel("Std Dev")
            fig1.savefig("StdDev_vs_TSYS_{0}.png".format(telescopes[0]))
            ax2.plot(tbl['{0}_TSYS'.format(telescopes[1])],
                     tbl['{0}_STDDEV'.format(telescopes[1])],
                     ',', alpha=0.8)
            ax2.set_xlabel("TSYS")
            ax2.set_ylabel("Std Dev")
            pl.draw()
            pl.show()
            fig2.savefig("StdDev_vs_TSYS_{0}.png".format(telescopes[1]))

    return datadict,tbldict
