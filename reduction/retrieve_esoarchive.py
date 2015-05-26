from astroquery.eso import Eso
import paths

# configure this
rawpath = '/scratch/aginsbur/apex/raw/'
username = 'aginsburg'

# don't touch the rest
Eso.cache_location = rawpath

Eso.login(username)

Eso.ROW_LIMIT = 1000000
tbl = Eso.query_instrument('apex', pi_coi='ginsburg', cache=False)
#stbl = tbl[np.char.startswith(tbl['Object'], 'Map') & (tbl['Scan Mode']=='OTF') & (tbl['Number of subscans'] > 10)]
#programs = set(tbl['ProgId'])
projids = set(tbl['APEX Project ID'])

for proj in projids:
    tbl = Eso.query_apex_quicklooks(proj[:-5])
    print(tbl)
    if proj[-4:] == '2015':
        files = Eso.retrieve_data(tbl['Product ID'])
