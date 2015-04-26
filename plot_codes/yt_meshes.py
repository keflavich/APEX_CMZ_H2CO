import numpy as np
import yt
import paths
from temperature_cubes import tcube,tcube_direct
from masked_cubes import cube303m
from spectral_cube import SpectralCube
from astropy.io import fits

tcube_direct._header['BTYPE'] = 'Temperature'
tcube_direct._header['BUNIT'] = 'K'
cube303m._header['BTYPE'] = 'flux'
#hdulist = fits.HDUList([cube303m.hdu, tcube_direct.hdu])
#hdu  = cube303m.hdu
#hdu.data = np.rollaxis(np.array([cube303m.hdu.data, tcube_direct.hdu.data]),0,3)


#ytds = yt.load(fits.HDUList(hdu), nan_mask=0.0, z_axis_decomp=True)
hdu1 = cube303m.hdu
hdu1.header['BTYPE'] = 'Brightness'
hdu1.header['BUNIT'] = 'K'
hdu2 = tcube_direct.hdu
hdu2.header['BTYPE'] = 'tkin'
hdu2.header['BUNIT'] = 'K'
ytds = yt.frontends.fits.FITSDataset(hdu1, auxiliary_files=[hdu2])

ytds.periodicity = (True,)*3

surface = ytds.surface(ytds.all_data(), ('fits','brightness'), 0.2)

surface.export_sketchfab(title='H2CO brightness colored by temperature',
                         description='none', color_field='tkin',
                         color_log=False, color_map='Blue-Red',
                         #color_field_min=20, color_field_max=300,
                        )

surface.export_obj('h2co_temperatures'.format(level),
                   color_field='tkin', color_log=False,
                   color_map='Blue-Red',
                   #color_field_min=20, color_field_max=300,
                   #transparency=0.5,
                  )

surface.export_sketchfab(title='H2CO brightness colored by temperature',
                         description='none', color_field='tkin',
                         color_log=False, color_map='Blue-Red',
                         #color_field_min=20, color_field_max=300,
                        )

for level in (0.1,0.2,0.3,0.4,0.5,0.6):


    surface = ytds.surface(ytds_sm.all_data(), ytds_sm.field_list[0], level)

    surface.export_ply('h2co_temperatures_cut{0:f}.ply'.format(level),
                       color_field='tkin', color_log=False,
                       color_map='Blue-Red',
                       #color_field_min=20, color_field_max=300,
                      )
