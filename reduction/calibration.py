import numpy as np
from astropy import log
from pyspeckit.spectrum.readers import read_class
from make_apex_cubes import all_apexfiles
import os
import paths
import collections

"""
Calibrators can be IDd with:
for ds in all_apexfiles:
     cl = read_class.ClassObject(ds)
     spdata = cl.get_spectra(line='CO(2-1)', range=np.array([-1,1,-1,1.])/206265.)
     rslt = set([h['OBJECT'] for d,h in spdata if 'SKY' not in h['OBJECT'] and 'COLD' not in h['OBJECT'] and 'TCAL' not in h['OBJECT'] and 'TSYS' not in h['OBJECT'] and 'TREC' not in h['OBJECT']])
     print os.path.basename(ds), rslt

H2CO cals:
for ds in all_apexfiles:
    cl = read_class.ClassObject(ds)
    spdata = cl.get_spectra(linere='shfi219|h2co', range=np.array([-1,1,-1,1.])/206265.)
    rslt = set([h['OBJECT'] for d,h in spdata if 'SKY' not in h['OBJECT'] and 'COLD' not in h['OBJECT'] and 'TCAL' not in h['OBJECT'] and 'TSYS' not in h['OBJECT'] and 'TREC' not in h['OBJECT'] and 'MAP' not in h['OBJECT']])
    print os.path.basename(ds), rslt
"""

calibrators = {'CO':('RAFGL1922', 'NGC6302', 'RT-SCO', 'IRAS15194-51',
                     'RAFGL2135', 'EP-AQR', 'PI1-GRU', 'RAFGL4211', 'CRL2688',
                     'IRAS17150-32'),
               'H2CO': ('SGRB2(M)', 'SGRB2(N)', ),
              }
all_tels = ['AP-H201-X202', 'AP-H201-X201', 'AP-H201-F102', 'AP-H201-F101']
telescopere = "|".join(all_tels)
linere = {'H2CO':'h2co|shfi219',
          'CO':'CO\(2-1\)'}

velocities = collections.defaultdict(lambda: 0.0, {'SGRB2(N)':65.0,
                                                   'SGRB2(M)':65.0,
                                                   'RAFGL2135':49.39,
                                                  })

cal_data = {species:{tel:{cal:[] for cal in calibrators[species]}
                     for tel in all_tels}
            for species in calibrators}
cal_spectra = {species:{tel:{cal:[] for cal in calibrators[species]}
                     for tel in all_tels}
            for species in calibrators}

for ds in all_apexfiles:
    cl = read_class.ClassObject(ds)
    sources = [s.strip() for s in cl.sources]

    for species in calibrators:
        for calibrator in calibrators[species]:
            if calibrator in sources:
                try:
                    spectra = cl.get_pyspeckit_spectra(source=calibrator,
                                                       linere=linere[species],
                                                       telescopere=telescopere,
                                                       range=np.array([-1,1,-1,1.])/206265.)
                except ValueError:
                    log.info("Skipped {0} {1}".format(calibrator, species))
                    continue
                vcen = float(spectra[0].header['VOFF'])
                if vcen == 0:
                    vcen = velocities[calibrator]
                vmin,vmax = vcen-50, vcen+50
                for sp in spectra:
                    if sp.data.max() > 1e4:
                        #log.info("Skipped {3} {1} {0} {2} - extreme values"
                        #         .format(sp.header['XTEL'], calibrator,
                        #                 sp.header['DOBS'], species))
                        continue
                    if species == 'H2CO': # use CH3OH because H2CO is absorbed
                        if 218.39258144781465 > sp.xarr.min() and 218.39258144781465 < sp.xarr.max():
                            sp.xarr.refX = 218.44005 #218.39258144781465 # CH3OH + 65.147 km/s
                        else:
                            sp.xarr.refX = 218.90336
                        sp.xarr.refX_units = 'GHz'
                    sp.xarr.convert_to_unit('km/s')
                    if sp.xarr.min() > 50 or sp.xarr.max() < -50:
                        log.info("Skipped {3} {1} {0} {2} - out of range"
                                 .format(sp.header['XTEL'], calibrator,
                                         sp.header['DOBS'], species))
                        continue
                    sp.specfit(guesses=[0,5,vcen,5], negamp=False,
                               fittype='vheightgaussian', xmin=vmin, xmax=vmax,
                               verbose=False)
                    while sp.specfit.parinfo.AMPLITUDE0.value == 0:
                        log.info("Repeating {3} {1} {0} {2} - zeros are unacceptable"
                                 .format(sp.header['XTEL'], calibrator,
                                         sp.header['DOBS'], species))
                        sp.specfit(guesses=[0,5,vcen,5], negamp=False,
                                   fittype='vheightgaussian', xmin=vmin, xmax=vmax,
                                   verbose=False)


                    rslt = (sp.header['DOBS'], sp.specfit.parinfo.AMPLITUDE0.value)
                    cal_data[species][sp.header['XTEL']][calibrator].append(rslt)
                    cal_spectra[species][sp.header['XTEL']][calibrator].append(sp)

                    log.info("{5} {1} {0} {2}: A={3}, v={4}"
                             .format(sp.header['XTEL'], calibrator,
                                     sp.header['DOBS'],
                                     sp.specfit.parinfo.AMPLITUDE0.value,
                                     sp.specfit.parinfo.SHIFT0.value,
                                     species))

                    if abs(sp.specfit.parinfo.SHIFT0.value - vcen) > 5:
                        log.warn("Velocity difference for {3} {0}:{1} was {2}"
                                 .format(sp.header['XTEL'], calibrator,
                                         abs(sp.specfit.parinfo.SHIFT0.value - vcen),
                                         species))
                
if not (os.path.exists(paths.fpath("calibration")) and
        os.path.isdir(paths.fpath('calibration'))):
    os.mkdir(paths.fpath('calibration'))

import pylab as pl
from astropy.time import Time

date1 = Time('2014-06-13 00:00:00.000', format='iso').jyear
date2 = Time('2014-04-23 00:00:00.000', format='iso').jyear

for species in cal_data:
    for tel in cal_data[species]:
        for source in cal_data[species][tel]:
            if (cal_data[species][tel][source]) == []:
                continue
            date,data = np.array(cal_data[species][tel][source]).T

            if len(data) > 0:

                pl.clf()
                pl.plot(date[data>0.5], data[data>0.5], '.')
                pl.xlabel("Decimal Year")
                pl.ylabel("Amplitude (K)")
                pl.title("{0} {1} {2}".format(species, tel, source))

                pl.savefig(paths.fpath('calibration/{0}_{1}_{2}_calibration.png'
                                       .format(species, tel, source)))

                mask = (date>date2) & (date<date1) & (data>1)
                mask2014 = (date>date1) & (data>1)
                if mask.sum() and mask2014.sum():
                    bad_mean = data[mask].mean()
                    good_mean = data[mask2014].mean()
                    print("{0:5s} {1:10s} {2:12s}: calfactor={3:5.3f} "
                          "npts(apr-june)={4:3d} npts(>june)={5:3d}".format(species,
                                                                     tel,
                                                                     source,
                                                                     good_mean/bad_mean,
                                                                     mask.sum(),
                                                                     mask2014.sum()))
