from astroquery.eso import Eso
import paths

def retrieve_ESO_files(rawpath='/scratch/aginsbur/apex/raw/',
                       username='aginsburg', projids=['O-085.F-9311A',
                                                      'E-085.B-0964A',
                                                      'E-093.C-0144A',
                                                      'E-095.C-0242A']):

    # don't touch the rest
    Eso.cache_location = rawpath

    Eso.login(username)

    Eso.ROW_LIMIT = 1000000
    #tbl = Eso.query_instrument('apex', pi_coi='ginsburg', cache=False)
    #stbl = tbl[np.char.startswith(tbl['Object'], 'Map') & (tbl['Scan Mode']=='OTF') & (tbl['Number of subscans'] > 10)]
    #programs = set(tbl['ProgId'])
    #projids = set(tbl['APEX Project ID'])

    for proj in projids:
        tbl = Eso.query_apex_quicklooks(proj)
        print(tbl)
        #files = Eso.retrieve_data(tbl['Product ID'])
