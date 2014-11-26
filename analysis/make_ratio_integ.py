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
    integ303321b = cube321.with_mask(bmasksm_rs).sum(axis=0)/cube303.with_mask(bmasksm_rs).sum(axis=0)
    hdu = cube321[0,:,:].hdu
    hdu.data = integ303321b.value
    hdu.writeto(paths.hpath('H2CO_321220_to_303202_integ_smoothmask.fits'))

                              
