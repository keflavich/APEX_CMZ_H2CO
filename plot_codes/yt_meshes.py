import yt
import paths



ytds_sm = yt.load(paths.hpath('APEX_H2CO_303_202_smooth_bl_snmasked.fits'),
               auxiliary_files=[paths.hpath('H2CO_321220_to_303202_cube_smooth_bl_temperature.fits')],
               nan_mask=0.0, z_axis_decomp=True)

ytds_sm.periodicity = (True,)*3

surface = ytds_sm.surface(ytds_sm.all_data(), ytds_sm.field_list[0], 0.2)

surface.export_sketchfab(title='H2CO brightness colored by temperature (smooth)',
                         description='none', color_field='tkin',
                         color_log=False, color_map='Blue-Red',
                         #color_field_min=20, color_field_max=300,
                        )

ytds = yt.load(paths.hpath('APEX_H2CO_303_202_bl_snmasked.fits'),
               auxiliary_files=[paths.hpath('H2CO_321220_to_303202_cube_bl_temperature.fits')],
               nan_mask=0.0, z_axis_decomp=True)

ytds.periodicity = (True,)*3

surface = ytds.surface(ytds.all_data(), ytds.field_list[0], 0.2)

surface.export_sketchfab(title='H2CO brightness colored by temperature',
                         description='none', color_field='tkin',
                         color_log=False, color_map='Blue-Red',
                         #color_field_min=20, color_field_max=300,
                        )

for level in (0.1,0.2,0.3,0.4,0.5,0.6):


    ytds_sm = yt.load(paths.hpath('APEX_H2CO_303_202_smooth_bl_snmasked.fits'),
                   auxiliary_files=[paths.hpath('H2CO_321220_to_303202_cube_smooth_bl_temperature.fits')],
                   nan_mask=0.0, z_axis_decomp=True)

    ytds_sm.periodicity = (True,)*3

    surface = ytds_sm.surface(ytds_sm.all_data(), ytds_sm.field_list[0], level)

    surface.export_ply('h2co_smooth_temperatures_cut{0:f}.ply'.format(level),
                       color_field='tkin', color_log=False,
                       color_map='Blue-Red',
                       #color_field_min=20, color_field_max=300,
                      )
