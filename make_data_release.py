"""
Copy data from the appropriate directories to the CDS directory for upload to A&A
"""
import os
from paths import h2copath,figurepath,hpath,rpath,fpath,mpath,molpath,tpath

integrated_files = [hpath(x) for x in ('APEX_H2CO_303_202_masked_moment0.fits',
                                       'APEX_H2CO_303_202_masked_smooth_moment0.fits',
                                       'APEX_H2CO_321_220_masked_moment0.fits',
                                       'APEX_H2CO_321_220_masked_smooth_moment0.fits',
                                       "H2CO_321220_to_303202_bl_integ_weighted.fits",
                                       "H2CO_321220_to_303202_bl_integ_masked_weighted.fits",
                                       "H2CO_321220_to_303202_bl_integ_temperature_dens1e4.fits",
                                       "H2CO_321220_to_303202_bl_integ_temperature_dens1e4_abund1e-10.fits",
                                       "H2CO_321220_to_303202_bl_integ_temperature_dens1e4_abund1e-8.fits",
                                       "H2CO_321220_to_303202_bl_integ_temperature_dens1e5.fits",
                                       "H2CO_321220_to_303202_bl_integ_temperature_dens3e4.fits",
                                       "H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e4.fits",
                                       "H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e4_abund1e-10.fits",
                                       "H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e4_abund1e-8.fits",
                                       "H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e5.fits",
                                       "H2CO_321220_to_303202_bl_integ_weighted_temperature_dens3e4.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_masked_weighted_temperature_dens1e4.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_masked_weighted_temperature_dens1e4_abund1e-10.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_masked_weighted_temperature_dens1e4_abund1e-8.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_masked_weighted_temperature_dens1e5.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_masked_weighted_temperature_dens3e4.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e4.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e4_abund1e-10.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e4_abund1e-8.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_weighted_temperature_dens1e5.fits",
                                       "pv_H2CO_321220_to_303202_bl_integ_weighted_temperature_dens3e4.fits",
                                      )
                   ]
cubes = [hpath(x)
         for x in ("H2CO_321220_to_303202_cube_bl.fits",
                   "H2CO_321220_to_303202_cube_smooth_bl.fits",
                   "APEX_H2CO_303_202_bl.fits",
                   "APEX_H2CO_321_220_bl.fits",
                   "APEX_H2CO_322_221_bl.fits",
                  )
        ] + [mpath(x)
             for x in ("APEX_13CO_2014_merge.fits",
                       "APEX_C18O_2014_merge.fits",
                       "APEX_H2CO_merge_high_plait_all.fits",
                      )
        ] + [molpath(x)
             for x in
             ("APEX_SiO_54_bl.fits",)
        ]


dendrograms = [hpath(x) for x in
               ("DendroMask_H2CO303202.hdf5",)
              ] + [tpath(x) for x in
                   ("fitted_line_parameters_Chi2Constraints.ipac",
                    "PPV_H2CO_Temperature.ipac",
                   )
              ]

if not os.path.isdir('cds'):
    os.mkdir('cds')
if not os.path.isdir('cds/integrated'):
    os.mkdir('cds/integrated')
if not os.path.isdir('cds/cubes'):
    os.mkdir('cds/cubes')
if not os.path.isdir('cds/catalogs'):
    os.mkdir('cds/catalogs')

for fn in integrated_files:
    fullfn = os.path.realpath(fn)
    shortfn = os.path.basename(fn)
    outpath = os.path.join('cds/integrated', shortfn)
    if not os.path.isfile(outpath):
        print(fn,outpath)
        os.link(fn, outpath)

for fn in cubes:
    fullfn = os.path.realpath(fn)
    shortfn = os.path.basename(fn)
    outpath = os.path.join('cds/cubes', shortfn)
    if not os.path.isfile(outpath):
        print(fn,outpath)
        os.link(fn, outpath)

for fn in dendrograms:
    fullfn = os.path.realpath(fn)
    shortfn = os.path.basename(fn)
    outpath = os.path.join('cds/catalogs', shortfn)
    if not os.path.isfile(outpath):
        print(fn,outpath)
        os.link(fn, outpath)
