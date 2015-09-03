from astropy.table import Table
import paths

to_remove = ['v_cen', 'r303321', 'er303321', '_idx']
to_rename = {'ratio321303': 'ratio303321',
             'eratio321303': 'eratio303321',
            }


flp = Table.read(paths.tpath('fitted_line_parameters_Chi2Constraints.ipac'), format='ascii.ipac')

for col in to_remove:
    if col in flp.colnames:
        print("removing {0} from flp".format(col))
        flp.remove_column(col)
for col in to_rename:
    if col in flp.colnames:
        print("Renaming {0} to {1} in flp".format(col, to_rename[col]))
        flp.rename_column(col, to_rename[col])


p2t = Table.read(paths.tpath('PPV_H2CO_Temperature_orbit.ipac'), format='ascii.ipac')

for col in to_remove:
    if col in p2t.colnames:
        print("removing {0} from p2t".format(col))
        p2t.remove_column(col)
for col in to_rename:
    if col in p2t.colnames:
        print("Renaming {0} to {1} in p2t".format(col, to_rename[col]))
        p2t.rename_column(col, to_rename[col])



tbl1 = tbl2 = False
cols1, cols2 = [],[]
with open(paths.tpath('README.rst'),'r') as f:
    for ii,line in enumerate(f.readlines()):
        if not line.strip():
            continue
        if 'fitted_line_parameters_Chi2Constraints_toCDS.ipac' in line:
            print("tbl1 starts: {0}".format(ii))
            tbl1=True
            tbl2=False
            continue
        elif 'PPV_H2CO_Temperature_orbit_toCDS.ipac' in line:
            print("tbl2 starts: {0}".format(ii))
            tbl2=True
            tbl1=False
            continue
        if tbl1:
            colname = line.split()[0]
            cols1.append(colname)
            if colname not in flp.colnames:
                print("{0} not in fitted_line_parameters_Chi2Constraints_toCDS.ipac".format(colname))
        elif tbl2:
            colname = line.split()[0]
            cols2.append(colname)
            if colname not in p2t.colnames:
                print("{0} not in PPV_H2CO_Temperature_orbit_toCDS.ipac".format(colname))

flpb = Table([flp[cn] for cn in cols1])
p2tb = Table([p2t[cn] for cn in cols2])

flpb.write(paths.tpath('fitted_line_parameters_Chi2Constraints_toCDS.ipac'), format='ascii.ipac')
p2tb.write(paths.tpath('PPV_H2CO_Temperature_orbit_toCDS.ipac'), format='ascii.ipac')
