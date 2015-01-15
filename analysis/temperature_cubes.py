import numpy as np
from paths import hpath, apath
from spectral_cube import SpectralCube, BooleanArrayMask
from ratio_cubes import (ratio303321, eratio303321, noise_flat, mask, ratioOK,
                         ratio303321sm, ratioOKsm, masksm)
from masked_cubes import (cube303m,cube321m,cube303msm,cube321msm,
                          cube303,cube321,cube303sm,cube321sm,
                          sncube, sncubesm)
from piecewise_rtotem import pwtem


tcubedata = np.empty(cube303m.shape)
tcubedata[~mask] = np.nan
tcubedata[mask] = pwtem(ratio303321[ratioOK])

tcube = SpectralCube(data=tcubedata, wcs=cube303m.wcs,
                     mask=cube303m.mask, meta={'unit':'K'},
                     header=cube303m.header,
                    )
assert tcube.header['CUNIT3'] == 'km s-1'
tcube_direct = tcube

# smoothed version
tcubedatasm = np.empty(cube303msm.shape)
tcubedatasm[~masksm] = np.nan
tcubedatasm[masksm] = pwtem(ratio303321sm[ratioOKsm])

tcubesm = SpectralCube(data=tcubedatasm, wcs=cube303msm.wcs,
                       mask=cube303msm.mask, meta={'unit':'K'},
                       header=cube303msm.header,)
assert tcubesm.header['CUNIT3'] == 'km s-1'
tcubesm_direct = tcubesm

for sm in ("","_smooth",'_321','_321smooth'):
    outpath = 'TemperatureCube_DendrogramObjects{0}_Piecewise.fits'.format(sm)
    tcube = SpectralCube.read(hpath(outpath))
    locals()['tcube_dend{0}'.format(sm)] = tcube
