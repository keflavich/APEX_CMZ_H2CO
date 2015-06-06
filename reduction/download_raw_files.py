from __future__ import print_function
import os
import retrieve_esoarchive
import cgi
import requests
from astropy.utils.data import download_files_in_parallel
from astropy.utils.console import ProgressBar

rawpath='/scratch/aginsbur/apex/raw/'
eso_username='aginsburg'

mpi_raw_files = ['https://dataverse.harvard.edu/api/v1/access/datafile/2509801',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2504422',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2509797',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2509803',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2509790',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2504424',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2504425',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2505538',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2509798',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2505539',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2509802',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2509804',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2509800',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2509799',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2505540',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2510176',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2505541',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2510177',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2505542',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2504384',
                 'https://dataverse.harvard.edu/api/v1/access/datafile/2504383']

def download_file(url, outdir=rawpath):
    r = requests.get(url, verify=False, stream=True)
    _, cdisp = cgi.parse_header(r.headers['content-disposition'])
    outfilename = cdisp['filename']
    fullname = os.path.join(outdir, outfilename)

    pb = ProgressBar(int(r.headers['content-length']))

    with open(fullname, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            f.write(chunk)
            f.flush()
            pb.update(pb._current_value + 1024)

    return fullname

def download_all():

    retrieve_esoarchive.retrieve_ESO_files(rawpath=rawpath, username=eso_username)

    for ii,fn in enumerate(mpi_raw_files):
        name = download_file(fn)
        print("Completed {2}: {0} of {1} files".format(ii, len(mpi_raw_files), name))
