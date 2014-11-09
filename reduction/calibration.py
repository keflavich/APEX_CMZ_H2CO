from pyspeckit.spectrum.readers import read_class
from make_apex_cubes import all_apexfiles

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
               'H2CO': ('SGRB2(M)', 'SGRB2(N)', 'SATURN', 'MARS'),
              }
all_tels = ['AP-H201-X202', 'AP-H201-X201', 'AP-H201-F102', 'AP-H201-F101']
linere = {'H2CO':'h2co|shfi219',
          'CO':'CO(2-1)'}

velo_info = {'PI1-GRU': [-50, 50, 0],
             'NGC6302': [-50, 50, 0],
             'RAFGL1922': [-50, 50, 0],
             'RAFGL2135': [-50, 50, 0],
             'EP-AQR': [-50, 50, 0],
             'IRAS15194-51': [-50, 50, 0],
            }

cal_data = {species:{tel:{cal:[] for cal in calibrators[species]}
                     for tel in all_tels}
            for species in calibrators}

for ds in all_apexfiles:
    cl = read_class.ClassObject(ds)

    for species in calibrators:
        for calibrator in calibrators[species]:
            if calibrator in cl.sources:
                spectra = cl.get_pyspeckit_spectra(source=calibrator,
                                                   linere=linere[species],
                                                   range=np.array([-1,1,-1,1.])/206265.)
                vmin,vmax,vcen = velo_info[calibrator]
                for sp in spectra:
                    sp.xarr.convert_to_unit('km/s')
                    sp.specfit(xmin=vmin, xmax=vmax)

                    rslt = (sp.header['DOBS'], sp.specfit.parinfo.AMPLITUDE0.value)
                    cal_data[species][sp.header['XTEL']][calibrator] = rslt

                    log.info("{4} {1} {0} {2}: {3}"
                             .format(sp.header['XTEL'], calibrator,
                                     sp.header['DOBS'],
                                     sp.specfit.parinfo.AMPLITUDE0.value,
                                     species))

                    if abs(sp.specfit.SHIFT0.value - vcen) > 5:
                        log.warn("Velocity difference for {3} {0}:{1} was {2}"
                                 .format(sp.header['XTEL'], calibrator,
                                         abs(sp.specfit.SHIFT0.value - vcen),
                                         species))
                

