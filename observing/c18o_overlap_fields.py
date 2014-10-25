"""
Find fields that were observed with the 219 GHz tuning, and find out if they
have been re-observed with the 218.9 GHz tuning.

This is mostly adequate, though there are unfortunately regions where even the
218.9 GHz tuning will result in some overlap.  It is not presently practical to
redo those.
"""
from astropy.table import Table
from paths import opath, rpath
import pyregion

tbl = Table.read(opath('observing_logs.ipac'), format='ipac')
co_good = co_ok = tbl['restf'] == 218900
co_bad = tbl['restf'] == 219000
maps_bad = [s.strip() for row in tbl[co_bad] for s in row['sources'].split(", ")]
maps_good = [s.strip() for row in tbl[co_good] for s in row['sources'].split(", ")]
bad_only = [m for m in maps_bad if m not in maps_good]

map_regions = pyregion.open(rpath('target_fields_8x8.reg'))

bad_regions = [m for m in map_regions if m.attr[1]['text'].split()[0].upper() in bad_only]
not_good_regions = [m for m in map_regions if m.attr[1]['text'].split()[0].upper() in maps_bad]
marked_regions = bad_regions + not_good_regions
for m in marked_regions:
    if m in bad_regions:
        m.attr[1]['color'] = 'red'
    else:
        m.attr[1]['color'] = 'green'
marked_regions = list(set(marked_regions))
pyregion.ShapeList(marked_regions).write(rpath("co_c18o_overlap.reg"))

for restf in [218.8, 218.9, 219.0]*u.GHz:
    print restf
    print "{0:18s} {1:20s}".format('VLSR(line)','V(12CO)_C18Oframe')
    for v in [-200,-100,0,100,200]*u.km/u.s:
        print "{0:18s} {1:20s}".format(v, -((((2 * (restf+6*u.GHz)) - 230.538*u.GHz*(1-v/c)) - 219.56036*u.GHz)/(219.56036*u.GHz) * c).to(u.km/u.s))
