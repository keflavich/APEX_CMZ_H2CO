from shfi_otf_pipeline import make_apex_cubes

datasets = {2013:make_apex_cubes.june2013datapath+'M-091.F-0019-2013-2013-06-11',
            'ao':make_apex_cubes.aorawpath+'E-085.B-0964A-2010_merge',
            2014:make_apex_cubes.april2014path+'E-093.C-0144A.2014APR03/E-093.C-0144A-2014-2014-04-02',
           }


frs = {}
all_headers = {}
for year,dataset in datasets.items():
    frs[year] = {}
    all_headers[year] = {}
    for lowhigh in ('low','high'):
        (spectra, headers, indices, data, hdrs,
         gal) = make_apex_cubes.load_dataset_for_debugging(skip_data=True,
                                                           lowhigh=lowhigh,
                                                           sourcename=None,
                                                           shapeselect=None,
                                                           backend='fft' if year=='ao' else 'xffts',
                                                           xscan=None,
                                                           dataset=dataset,
                                                           datapath='')
        frs[year][lowhigh] = make_apex_cubes.hdr_to_freq(hdrs[0])
        all_headers[year][lowhigh] = headers

print {(x,y):(round(frs[x][y].min()/1e3,1),round(frs[x][y].max()/1e3,1)) for x in frs for y in frs[x]}
"""
{(2013, 'high'): (217.5, 220.0),
 (2013, 'low'):  (216.0, 218.5),
 (2014, 'high'): (218.4, 220.9),
 (2014, 'low'):  (216.9, 219.4),
 ('ao', 'high'): (218.0, 219.0),
 ('ao', 'low'):  (216.8, 217.8)}
"""
