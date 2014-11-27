import numpy as np
import paths
from ratio_cubes import (ratio303321, eratio303321,
                         ratio303321sm, eratio303321sm,
                         ratiocube_303321,
                         ratiocubesm_303321,
                        )
from masked_cubes import (cube303, cube303sm, cube303m, cube321m, cube303msm,
                          cube321msm, cube321, cube321sm, bmask, bmasksm, mask,
                          masksm, bmasksm_rs)
from noise import (noise, noise_cube, sm_noise, cube303nm, cube303nmsm,
                   cube321nm, cube321nmsm)


def make_ratio_integ():
    npixflat = bmasksm_rs.include().sum(axis=0)
    sum321 = cube321.with_mask(bmasksm_rs).sum(axis=0)
    sum303 = cube303.with_mask(bmasksm_rs).sum(axis=0)
    integ303321b = sum321/sum303
    error303 = cube303.with_mask(~bmasksm_rs).std(axis=0)
    error321 = cube321.with_mask(~bmasksm_rs).std(axis=0)

    eintegratio303321 = ((integ303321b.value**2 * (error303.value**2/(sum303.value/npixflat)**2 +
                                                   error321.value**2/(sum321.value/npixflat)**2))**0.5)

    ok = (integ303321b.value/eintegratio303321) > 3
    integ303321ok = integ303321b.value
    integ303321ok[~ok] = np.nan

    hdu = cube321[0,:,:].hdu
    hdu.data = integ303321ok
    hdu.writeto(paths.hpath('H2CO_321220_to_303202_integ_smoothmask.fits'),clobber=True)
    hdu.data = eintegratio303321
    hdu.writeto(paths.hpath('H2CO_321220_to_303202_einteg_smoothmask.fits'),clobber=True)


def make_ratio_max():
    npixflat = bmasksm_rs.include().sum(axis=0)
    max321 = cube321.with_mask(bmasksm_rs).max(axis=0)
    max303 = cube303.with_mask(bmasksm_rs).max(axis=0)
    maxratio303321b = max321/max303
    error303 = cube303.with_mask(~bmasksm_rs).std(axis=0)
    error321 = cube321.with_mask(~bmasksm_rs).std(axis=0)

    eintegratio303321 = ((maxratio303321b.value**2 * (error303.value**2/(max303.value)**2 +
                                                      error321.value**2/(max321.value)**2))**0.5)

    ok = (maxratio303321b.value/eintegratio303321) > 3
    maxratio303321ok = maxratio303321b.value
    maxratio303321ok[~ok] = np.nan

    hdu = cube321[0,:,:].hdu
    hdu.data = maxratio303321ok
    hdu.writeto(paths.hpath('H2CO_321220_to_303202_maxratio_smoothmask.fits'),clobber=True)
    hdu.data = eintegratio303321
    hdu.writeto(paths.hpath('H2CO_321220_to_303202_emaxratio_smoothmask.fits'),clobber=True)


