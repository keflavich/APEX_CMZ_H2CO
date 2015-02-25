from astropy import table
from astropy import log
import numpy as np

def flag_dendro(dend, catalog=None, smooth=False, pixels_with_bad=[],
                colname='IsNotH2CO', flag_descendants=True):
    """
    Remove / flag as "bad" objects in the dendrogram that are not H2CO

    This was done manually; there is now (Jan 17, 2015) an automated approach
    in dendro_temperature
    """

    if catalog is not None and colname not in catalog.colnames:
        catalog.add_column(table.Column(name=colname, dtype=bool,
                                        data=np.zeros(len(catalog),
                                                      dtype='bool')))
    for obj in dend:
        obj.bad = False

    for x,y,z in pixels_with_bad:
        z = z-64 # Jan 25: the cubes were shifted
        if z < 0: continue
        bad_obj = dend.structure_at([z/(2 if smooth else 1),y,x])
        if bad_obj and flag_descendants:
            bad_obj = bad_obj.ancestor
            bad_obj.bad = True
            if catalog is not None:
                catalog[bad_obj.idx][colname] = True
            for obj in bad_obj.descendants:
                obj.bad = True
                catalog[obj.idx][colname] = True
        elif bad_obj:
            obj = bad_obj
            # flag the ancestors, but not the ancenstors' parents
            while obj:
                obj.bad = True
                if catalog is not None:
                    catalog[obj.idx][colname] = True
                obj = obj.parent
        else:
            # sometimes these don't get IDd?
            pass

def flag_hc3n(dend, catalog, smooth):
    pixels_with_bad = [(518,122,142),
                       (508,123,10),
                       (533,124,48),
                       (857,103,126),
                       (904,102,105),
                       (884,95,118),
                       (515, 108, 1),
                      ]
    flag_dendro(dend, catalog, smooth, pixels_with_bad=pixels_with_bad,
                colname='IsNotH2CO')

    if issubclass(catalog['IsNotH2CO'].dtype.type, str):
        col = catalog['IsNotH2CO']
        catalog.remove_column('IsNotH2CO')
        catalog.add_column(table.Column(name='IsNotH2CO',
                                        dtype='bool',
                                        data=col=='True'))


    log.info("Flagged {0} dendrogram objects as HC3N".format(catalog['IsNotH2CO'].sum()))

def flag_absorption(dend, catalog=None, smooth=False):
    """
    Flag out things associated with absorption in one or both of the Sgr B2
    lines
    """
    # add 64 because I got these from the post-Jan25 cubes
    pixels_with_bad = [(511,124,221+64),
                       (512,124,181+64),
                       (516,125,218+64),
                       (517,121,197+64),
                       (515,120,240+64),
                      ]
    flag_dendro(dend, catalog, smooth, pixels_with_bad=pixels_with_bad,
                colname='IsAbsorption', flag_descendants=False)

    if issubclass(catalog['IsAbsorption'].dtype.type, str):
        col = catalog['IsAbsorption']
        catalog.remove_column('IsAbsorption')
        catalog.add_column(table.Column(name='IsAbsorption',
                                        dtype='bool',
                                        data=col=='True'))

    log.info("Flagged {0} dendrogram objects as Absorption".format(catalog['IsAbsorption'].sum()))
