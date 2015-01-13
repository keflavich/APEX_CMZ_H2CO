from astroquery.eso import Eso
import numpy as np
Eso.ROW_LIMIT = 1000000
tbl = Eso.query_instrument('apex', pi_coi='ginsburg', cache=False)
stbl = tbl[np.char.startswith(tbl['Object'], 'Map') & (tbl['Scan Mode']=='OTF') & (tbl['Number of subscans'] > 10)]

for object in np.unique(stbl['Object']):
    matches = stbl['Object'] == object
    total_exptime = stbl['ExpTime'][matches].sum()
    files = np.unique(stbl['ORIGFILE'][matches])
    print object,"\t",total_exptime,"\t","\t".join([x[11:21] for x in files])
    
dates = [x[11:21] for x in stbl['ORIGFILE']]

for d in np.unique(dates):
    matches = np.array(dates)==d
    objects = np.unique(stbl[matches]['Object'])
    print d, tuple(list(objects))
