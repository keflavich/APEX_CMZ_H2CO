"""
Script to convert ds9 region files to APEX catalog files
"""
import pyregion
import numpy as np
from astropy import units as u
from astropy import coordinates
import time

if False:
    print "APEX: Old... these were the original pointings from the June 2013 run",
    t0=time.time()
    with open('target_fields_reorganized.APEX','w') as outf:

        r = pyregion.open('target_fields_reorganized.reg')
        for ii,reg in enumerate(r):
            # only observe new target fields
            if reg.attr[1]['color'] == 'orange':
                c = coordinates.ICRS(reg.coord_list[0],reg.coord_list[1],unit=(u.deg,u.deg))
                print >>outf,"Map_{:03}   EQ 2000  {:15} {:15}  LSR 0.0  ! comment".format(ii+1,c.ra.to_string(unit=u.h, sep=':'),c.dec.to_string(unit=u.deg,sep=':'))
    print "Done (%0.1g s)" % (time.time()-t0)

    print "APEX: P93, PA = 0, 4x4... these are NOT USED but were entertained as an early option",
    t0=time.time()
    with open('target_fields_P93ESO.APEX','w') as outf:

        r = pyregion.open('observed_fields_tbobs_fk5.reg')
        for ii,reg in enumerate(r):
            # only observe new target fields
            if reg.attr[1]['color'] == 'orange':
                c = coordinates.ICRS(reg.coord_list[0],reg.coord_list[1],unit=(u.deg,u.deg))
                print >>outf,"P93_Map_{:03}   EQ 2000  {:15} {:15}  LSR 0.0  ! comment".format(ii+1,c.ra.to_string(unit=u.h, sep=':'),c.dec.to_string(unit=u.deg,sep=':'))
    print "Done (%0.1g s)" % (time.time()-t0)


print "APEX: P93, PA = 58, 8x8... these are the 2014 targets!  They are 8x8' grids aligned with the Galactic plane",
print "The numbering scheme is somewhat arbitrary, and is actually due to an error in the logic below."
print "0-32 are 'top priority'"
print "55-58 are 'intermediate priority'"
print "115-124 are 'tertiary priority'"
print "Total fields = 42."
t0=time.time()
with open('target_fields_P93ESO_8x8.APEX','w') as outf:
    apex_reg_template = "Map_{:03}   EQ 2000  {:15} {:15}  LSR 0.0  ! {}"

    r = pyregion.open('target_fields_8x8.reg')
    off_regions = coordinates.ICRS(*np.array(zip(*[(x.coord_list[0],x.coord_list[1])
                                             for x in
                                             pyregion.open('off_positions_selectedfromDame2001.reg')
                                             if x.attr[1]['color'] == 'green']))*u.deg)
    #print
    for ii,reg in enumerate(r):
        # only observe new target fields
        #print ii+1,reg,reg.attr[1]['color']
        if reg.attr[1]['color'] == 'blue':
            c = coordinates.ICRS(reg.coord_list[0],reg.coord_list[1],unit=(u.deg,u.deg))
            closest_off = np.argmin(off_regions.separation(c))
            print >>outf, apex_reg_template.format(ii+1,
                                                   c.ra.to_string(unit=u.h,
                                                                  sep=':'),
                                                   c.dec.to_string(unit=u.deg,
                                                                   sep=':'),
                                                   "Top Priority.  Off %i" % (closest_off+1))

    offsetnum = ii


    for ii,reg in enumerate(r):
        # Secondary Priority
        if reg.attr[1]['color'] == 'cyan':
            c = coordinates.ICRS(reg.coord_list[0],
                                 reg.coord_list[1],unit=(u.deg,u.deg))
            closest_off = np.argmin(off_regions.separation(c))
            print >>outf, apex_reg_template.format(offsetnum+ii+1,
                                                   c.ra.to_string(unit=u.h,
                                                                  sep=':'),
                                                   c.dec.to_string(unit=u.deg,
                                                                   sep=':'),
                                                   "Secondary Priority.  Off %i" % (closest_off+1))

    offsetnum += ii

    for ii,reg in enumerate(r):
        # Tertiary priority
        if reg.attr[1]['color'] == 'purple':
            c = coordinates.ICRS(reg.coord_list[0],
                                 reg.coord_list[1],unit=(u.deg,u.deg))
            closest_off = np.argmin(off_regions.separation(c))
            print >>outf, apex_reg_template.format(offsetnum+ii+1,
                                                   c.ra.to_string(unit=u.h,
                                                                  sep=':'),
                                                   c.dec.to_string(unit=u.deg,
                                                                   sep=':'),
                                                   "Tertiary Priority.  Off %i" % (closest_off+1))
print "Done (%0.1g s)" % (time.time()-t0)
