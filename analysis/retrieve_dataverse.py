import requests
import os
import sys

#dv = requests.get('https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/RYEANM')
dv = requests.get('https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/27601')
dv.raise_for_status()

data = dv.json()

files = data['data']['latestVersion']['files']

root = '/Volumes/external/apex/'

for fd in files:

    result = requests.get('https://dataverse.harvard.edu/api/access/datafile/{id}'
                          .format(**fd['dataFile']),
                          stream=True)
    result.raise_for_status()

    with open(os.path.join(root, fd['dataFile']['filename']), 'wb') as fh:
        for chunk in result.iter_content(chunk_size=1024**2):
            fh.write(chunk)
            sys.stdout.write('.')
