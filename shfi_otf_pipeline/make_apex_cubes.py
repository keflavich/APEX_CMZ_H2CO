import pyspeckit
from pyspeckit.spectrum.readers import read_class
from pyspeckit import cubes
import numpy as np
from astropy import wcs
from astropy import coordinates
from astropy import units as u
from astropy import constants
from .progressbar import ProgressBar
from astropy.convolution import convolve, Gaussian1DKernel
from sdpy import makecube
from astropy.io import fits
from FITS_tools import cube_regrid
from FITS_tools.load_header import get_cd
from astropy.wcs import WCS
import FITS_tools
import scipy.ndimage
import scipy.linalg
import time
import mpl_plot_templates
import pylab as pl
import os
import errno
from astropy import log
import glob
from scipy.ndimage import filters
from scipy import signal,interpolate
import warnings
import image_tools
import spectral_cube

datasets_ao = ['O-085.F-9311A-2010','E-085.B-0964A-2010']
datasets_2013 = ['M-091.F-0019-2013-2013-06-08',
                 'M-091.F-0019-2013-2013-06-11',
                 'M-091.F-0019-2013-2013-06-12',
                 'M-091.F-0019-2013-2013-06-13']
datasets_2014 = {'E-093.C-0144A.2014JUN01/E-093.C-0144A-2014-2014-05-31': ('MAP_007',),
                 'E-093.C-0144A.2014MAY30/E-093.C-0144A-2014-2014-05-29': ('MAP_002','MAP_003','MAP_004'),
                 'E-093.C-0144A.2014MAY31/E-093.C-0144A-2014-2014-05-30': ('MAP_005','MAP_006'),
                 'E-093.C-0144A.2014APR02/E-093.C-0144A-2014-2014-04-01': ('MAP_001',),
                 'E-093.C-0144A.2014APR03/E-093.C-0144A-2014-2014-04-02': ('MAP_001',),
                 'E-093.C-0144A.2014JUN02/E-093.C-0144A-2014-2014-06-01': ('MAP_009','MAP_010','MAP_008',),
                 'E-093.C-0144A.2014JUN03/E-093.C-0144A-2014-2014-06-02': ('MAP_011','MAP_012','MAP_013', 'MAP_018', 'MAP_019'),
                 'E-093.C-0144A.2014JUN06/E-093.C-0144A-2014-2014-06-05': ('Map_020', 'Map_021', 'Map_022', 'Map_023', 'Map_024', 'Map_025'),
                 # There is some corrupt data in 06-06
                 'E-093.C-0144A.2014JUN07/E-093.C-0144A-2014-2014-06-06': ('Map_001', 'Map_026', 'Map_027', 'Map_028', 'Map_029', 'Map_030'),
                 'E-093.C-0144A.2014JUL29/E-093.C-0144A-2014-2014-07-28': ('MAP_002', 'MAP_001'),
                 'E-093.C-0144A.2014JUL29/E-093.C-0144A-2014-2014-07-29': ('MAP_002',),
                 'E-093.C-0144A.2014JUL30/E-093.C-0144A-2014-2014-07-29': ('MAP_004', 'MAP_002', 'MAP_003'),
                 'E-093.C-0144A.2014JUL31/E-093.C-0144A-2014-2014-07-30': ('MAP_005', 'MAP_006'),
                 'E-093.C-0144A.2014AUG01/E-093.C-0144A-2014-2014-07-31': ('MAP_006', 'MAP_007', 'MAP_008', 'MAP_009', 'MAP_012', 'MAP_011', 'MAP_010'),
                 'E-093.C-0144A.2014AUG01/E-093.C-0144A-2014-2014-08-01': ('MAP_013',),
                 'E-093.C-0144A.2014AUG02/E-093.C-0144A-2014-2014-08-01': ('MAP_013', 'MAP_018'),
                 'E-093.C-0144A.2014AUG09/E-093.C-0144A-2014-2014-08-07': ('MAP_024', 'MAP_022', 'MAP_023', 'MAP_025'),
                 'E-093.C-0144A.2014AUG09/E-093.C-0144A-2014-2014-08-08': ('MAP_027', 'MAP_026'),
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-04-24': ('MAP_115','MAP_116',),
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-04-30': ('MAP_116',),
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-02': ('MAP_116',),
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-04': ('MAP_115','MAP_116',),
                 # should be 05-07: map117
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-08': ('MAP_117','MAP_118',),
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-09': ('MAP_119','MAP_118',),
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-10': ('MAP_120','MAP_121','MAP_119',),
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-11': ('MAP_121','MAP_122','MAP_123','MAP_124',),
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-12': ('MAP_055','MAP_056','MAP_124',),
                 'M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-13': ('MAP_031','MAP_032','MAP_057','MAP_058',),
                 'M-093.F-0009-2014-2014-07-20':('MAP_024'),
                 'M-093.F-0009-2014-2014-07-19':['MAP_024'],
                 'M-093.F-0009-2014-2014-07-14':('MAP_024','MAP_025'),
                 'M-093.F-0009-2014-2014-07-13':['MAP_028', 'MAP_026', 'MAP_027', 'MAP_024', 'MAP_025'],
                 'M-093.F-0009-2014-2014-07-12':['MAP_028', 'MAP_029'],
                 'M-093.F-0009-2014-2014-07-11':['MAP_029', 'MAP_030'],
                 'M-093.F-0009-2014-2014-07-10':['MAP_031', 'MAP_030'],
                }
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-04-24 ['MAP_115', 'MAP_116']
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-04-30 ['MAP_116']
#E-093.C-0144A.2014JUN01/E-093.C-0144A-2014-2014-05-31 ['MAP_007']
#E-093.C-0144A.2014APR03/E-093.C-0144A-2014-2014-04-02 ['MAP_001']
#E-093.C-0144A.2014JUN02/E-093.C-0144A-2014-2014-06-01 ['MAP_008', 'MAP_009', 'MAP_010']
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-09 ['MAP_118', 'MAP_119']
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-08 ['MAP_118', 'MAP_117']
#E-093.C-0144A.2014JUN03/E-093.C-0144A-2014-2014-06-02 ['MAP_013', 'MAP_012', 'MAP_011', 'MAP_019', 'MAP_018']
#E-093.C-0144A.2014MAY31/E-093.C-0144A-2014-2014-05-30 ['MAP_005', 'MAP_006']
#E-093.C-0144A.2014MAY30/E-093.C-0144A-2014-2014-05-29 ['MAP_004', 'MAP_002', 'MAP_003']
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-02 ['MAP_116']
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-04 ['MAP_115', 'MAP_116']
#E-093.C-0144A.2014APR02/E-093.C-0144A-2014-2014-04-01 ['MAP_001']
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-10 ['MAP_119', 'MAP_121', 'MAP_120']
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-11 ['MAP_123', 'MAP_122', 'MAP_121', 'MAP_124']
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-12 ['MAP_124', 'MAP_056', 'MAP_055']
#M-093.F-0009-2014-2014-04/M-093.F-0009-2014-2014-05-13 ['MAP_031', 'MAP_058', 'MAP_057', 'MAP_032']


june2013datapath = '/Users/adam/work/h2co/apex/june2013/raw/M-091.F-0019-2013/'
june2013path = '/Users/adam/work/h2co/apex/june2013/'
april2014path = '/Users/adam/work/h2co/apex/april2014/'
h2copath = '/Users/adam/work/h2co/apex/h2co_cubes/'
mergepath = '/Users/adam/work/h2co/apex/merged_datasets/'
aorawpath = '/Users/adam/work/h2co/apex/2010_reduced/2010_raw/'
aopath = '/Users/adam/work/h2co/apex/2010_reduced/'
diagplotdir = '/Users/adam/work/h2co/apex/diagnostic_plots/'


def mkdir_p(path):
    """ http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def checkdir_makedir(path):
    dpath = os.path.split(path)[0]
    if not os.path.exists(dpath) and dpath:
        mkdir_p(dpath)


def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default

    """

    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m

def debug_and_load(test='test'):

    spectra,headers,indices,data,hdrs,gal = load_dataset_for_debugging(skip_data=False, lowhigh='high')

    make_blanks_freq(gal, hdrs[0], test, clobber=True)
    dmeansub,gal,hdrs = process_data(data, gal, hdrs, dataset=test,
                                     subspectralmeans=True, scanblsub=False)
    add_apex_data(dmeansub, hdrs, gal, test, retfreq=True, varweight=True,)
    dscube = cube_regrid.downsample_cube(fits.open(test+".fits")[0], factor=4)
    dscube.writeto(test+"_ds.fits",clobber=True)

    make_blanks_freq(gal, hdrs[0], test+"_blsub", clobber=True)
    dspecsub,gal,hdrs = process_data(data, gal, hdrs, dataset=test+"_blsub",
                                     subspectralmeans=True, scanblsub=True)
    add_apex_data(dspecsub, hdrs, gal, test+"_blsub", retfreq=True, varweight=True,)
    dscube = cube_regrid.downsample_cube(fits.open(test+"_blsub.fits")[0], factor=4)
    dscube.writeto(test+"_blsub_ds.fits",clobber=True)

    make_blanks_freq(gal, hdrs[0], test+"_pcasub", clobber=True)
    dpcasub,gal,hdrs = process_data(data, gal, hdrs, dataset=test+"_pcasub",
                                    subspectralmeans=True, scanblsub=True,
                                    pca_clean=True, pcakwargs={})
    add_apex_data(dpcasub, hdrs, gal, test+"_pcasub", retfreq=True, varweight=True,)
    dscube = cube_regrid.downsample_cube(fits.open(test+"_pcasub.fits")[0], factor=4)
    dscube.writeto(test+"_pcasub_ds.fits",clobber=True)

    freq = hdr_to_freq(hdrs[0])
    mask = make_line_mask(freq)

    return spectra,headers,indices,data,hdrs,gal,dspecsub,dmeansub,dpcasub,freq,mask

def load_dataset_for_debugging(lowhigh='low', downsample_factor=8,
                               dataset='M-091.F-0019-2013-2013-06-11',
                               datapath=june2013datapath,
                               xscan=37986,
                               sourcename='SGRA',
                               shapeselect=4096,
                               backend='xffts',
                               skip_data=True):
    """
    Example:

    spectra,headers,indices, data,hdrs,gal = load_dataset_for_debugging(skip_data=False)
    make_blanks_freq(gal, hdrs[0], 'test', clobber=True)
    noise = np.std(data,axis=1)
    freq_step = np.array([h['FRES'] for h in hdrs])
    exptime = np.array([h['EXPOSURE'] for h in hdrs])
    tsys = np.array([h['TSYS'] for h in hdrs])
    diagplot(data, tsys, noise, 'test')
    add_apex_data(data, hdrs, gal, cubefilename, retfreq=True, varweight=True,)
    """

    if lowhigh not in ('low','high'):
        raise ValueError
    if backend == 'xffts': 
        xtel = 'AP-H201-X202' if lowhigh=='low' else 'AP-H201-X201'
    else:
        xtel = 'AP-H201-F101' if lowhigh == 'high' else 'AP-H201-F102'

    apex_filename=datapath+dataset+".apex"

    spectra,headers,indices = load_apex_cube(apex_filename,
                                             skip_data=skip_data,
                                             downsample_factor=downsample_factor)
    data, hdrs, gal = select_apex_data(spectra, headers, indices,
                                       sourcename=sourcename,
                                       shapeselect=shapeselect,
                                       tsysrange=[100,325],
                                       xtel=xtel,
                                       rchanrange=None,
                                       xscan=xscan,
                                       skip_data=skip_data)

    return spectra,headers,indices, data,hdrs,gal

def get_sourcenames(headers):
    return list(set([h['SOURC'].strip() for h in headers]))

def load_apex_cube(apex_filename='data/E-085.B-0964A-2010_merge.apex',
                   skip_data=False, DEBUG=False, downsample_factor=None):
    spectra,headers,indices = read_class.read_class(apex_filename, start=1024,
                                                    DEBUG=DEBUG,
                                                    skip_data=skip_data,
                                                    downsample_factor=downsample_factor)

    for h,i in zip(headers,indices):
        h.update(i)

    return spectra,headers,indices

def add_apex_cube(apex_filename='data/E-085.B-0964A-2010_merge.apex',
                  cubefilename='APEX_H2CO_Ao', clobber=True,
                  kernel_fwhm=10./3600., **kwargs):
    spectra,headers,indices = load_apex_cube(apex_filename)

    data, hdrs, gal = select_apex_data(spectra, headers, indices, **kwargs)

    add_apex_data(data, hdrs, gal, cubefilename, kernel_fwhm=kernel_fwhm)

def select_apex_data(spectra,headers,indices, sourcename=None,
                     shapeselect=None, tsysrange=None, rchanrange=None,
                     xscan=None,
                     xtel=None,
                     skip_data=False,
                     galactic_coordinate_range=[[-2,2],[-2,2]]):

    log.info("Determining RA/Dec")
    ra,dec = zip(*[(h['RA']+h['RAoff']/np.cos(h['DEC']/180.*np.pi),
                    h['DEC']+h['DECoff']) for h in headers])
    log.info("Determining Galactic coordinates")
    gal = coordinates.SkyCoord(np.array(ra)*u.deg,
                               np.array(dec)*u.deg,
                               frame='icrs').galactic
    #gal.l.wrap_angle = 180*u.deg
    if galactic_coordinate_range is not None:
        (lmin,lmax),(bmin,bmax) = galactic_coordinate_range
        galOK = ((gal.l.wrap_at(180*u.deg).deg > lmin) &
                 (gal.l.wrap_at(180*u.deg).deg < lmax) &
                 (gal.b.deg > bmin) &
                 (gal.b.deg < bmax))
    else:
        galOK = True

    
    if isinstance(sourcename, (list,tuple)):
        sourceOK = np.array([h['SOURC'].strip() in sourcename for h in headers])
    elif sourcename is not None:
        sourceOK = np.array([h['SOURC'].strip()==sourcename for h in headers])
    else:
        sourceOK = True

    if xscan is not None:
        xscanOK = np.array([h['XSCAN']==xscan for h in headers])
    else:
        xscanOK = True


    if xtel is not None:
        xtelOK = np.array([h['XTEL'].strip()==xtel for h in headers])
    else:
        xtelOK = True


    if tsysrange is not None:
        tsys = np.array([h['TSYS'] for h in headers])
        tsysOK = (tsys>tsysrange[0]) & (tsys<tsysrange[1])
    else:
        tsysOK = True

    if rchanrange is not None:
        rchan = np.array([h['RCHAN'] if 'RCHAN' in h else np.inf for h in headers])
        rchanOK = (rchan>rchanrange[0]) & (rchan<rchanrange[1])
    else:
        rchanOK = True

    mostOK = galOK & sourceOK & tsysOK & rchanOK & xtelOK & xscanOK

    if not skip_data:
        log.info("Shaping data")
        data1 = np.array(spectra)
        shapes = np.array([d.shape for d in data1])
        if shapeselect is not None:
            OKshapes = (shapes == shapeselect).squeeze()
        elif len(np.unique(shapes[mostOK])) > 1:
            raise ValueError("Inconsistent shapes.")
        else:
            OKshapes = True
    else:
        OKshapes = True


    allOK = mostOK & OKshapes
    if allOK.sum() == 0:
        raise ValueError("Data selection yielded empty.  Sourcename={0}".format(sourcename))

    if skip_data:
        data = None
    else:
        data = np.array(data1[allOK].tolist())

    hdrs = [h for h,K in zip(headers,allOK) if K]
    gal = gal[allOK]

    return data,hdrs,gal

def process_data(data, gal, hdrs, dataset, scanblsub=False,
                 subspectralmeans=True, verbose=False, noisefactor=1.5,
                 linemask=False, automask=2,
                 zero_edge_pixels=0,
                 subtract_time_average=False,
                 pca_clean=False,
                 timewise_pca=True,
                 pcakwargs={},
                 **kwargs):

    timeaxis = 0
    freqaxis = 1

    log.info("Processing {0}".format(dataset))

    if zero_edge_pixels:
        # Force the Nth first/last frequency pixels to zero
        data[:,:zero_edge_pixels] = 0
        data[:,-zero_edge_pixels:] = 0

    # flag extremely bad pixels (don't know where these come from, scary!)
    extremely_bad = (data > 1e10) | (data < -1e10)
    # Set to zero rather than nan to avoid masking-related issues below
    data[extremely_bad] = 0

    if subspectralmeans:
        data = data - data.mean(axis=freqaxis)[:,None]

    obsids = np.array([h['XSCAN'] for h in hdrs])

    # for plotting and masking, determine frequency array
    freq = hdr_to_freq(hdrs[0])

    scans = identify_scans_fromcoords(gal)

    if scanblsub:

        data_diagplot(data, dataset+"_presub", scans=scans, freq=freq,
                      **kwargs)
        for ii,xscan in enumerate(np.unique(obsids)):
            match = obsids == xscan
            # maybe mask=mask_pix.max(axis=timeaxis), ?
            #mask=mask_pix[ii], 
            data_diagplot(data[match], dataset+"_presub_obs%i" % xscan,
                          freq=freq, **kwargs)

        if linemask:
            mask = make_line_mask(freq)
        else:
            mask = None
        dsub,mask_pix = subtract_scan_linear_fit(data, scans, mask_pixels=mask,
                                                 verbose=verbose,
                                                 automask=automask,
                                                 smooth_all=True,
                                                 return_mask=True)
        if len(mask_pix) == 0:
            mask = None
        else:
            mask = mask_pix.max(axis=timeaxis).astype('bool')
    elif subtract_time_average:
        # subtracting mean spectrum from all spectra
        dsub = data - data.mean(axis=timeaxis)
        mask = None
    else:
        mask = None
        dsub = data

    if pca_clean:
        t0 = time.time()
        if timewise_pca:
            dsub = PCA_clean(dsub.T, smoothing_scale=False,
                             diagplotfilename=os.path.join(diagplotdir,
                                                           dataset+"_time_pca_diagnostic.png"),
                             **pcakwargs).T
        else:
            # DON'T remove the mean: that's dealt with in 'spectral baselining' in
            # a more conservative fashion
            dmean = dsub.mean(axis=0)
            dsub = PCA_clean(dsub-dmean, 
                             diagplotfilename=os.path.join(diagplotdir,
                                                           dataset+"_pca_diagnostic.png"),
                             **pcakwargs) + dmean
        log.info("PCA cleaning took {0} seconds".format(time.time()-t0))

    # Standard Deviation can be fooled by obscene outliers
    #noise = MAD(dsub,axis=freqaxis)
    noise = np.std(dsub,axis=freqaxis)
    freq_step = np.array([h['FRES'] for h in hdrs])
    exptime = np.array([h['EXPOSURE'] for h in hdrs])
    tsys = np.array([h['TSYS'] for h in hdrs])
    theoretical_rms = 2.0**0.5*tsys/(np.abs(freq_step*1.0e6)*exptime)**0.5
    # extra factor 1.5 to avoid overflagging... don't know why really
    bad = noise > (theoretical_rms*noisefactor)

    # SgrB2 has higher noise.  Don't flag it out.
    sgrb2 = ((gal.l.wrap_at(180*u.deg).deg > 0.64) &
             (gal.l.wrap_at(180*u.deg).deg<0.68) &
             (gal.b.deg>-0.05) &
             (gal.b.deg<-0.001))
    bad[sgrb2] = False

    # pre-flagging diagnostic
    diagplot(dsub, tsys, noise, dataset+"_preflag", freq=freq, mask=mask,
             scans=scans, **kwargs)

    if np.count_nonzero(bad) == bad.size:
        import ipdb; ipdb.set_trace()
        raise ValueError("All data will be flagged out; something is amiss.")

    dsub = dsub[True-bad]
    obsids = obsids[True-bad]
    tsys = tsys[True-bad]
    noise = noise[True-bad]

    gal = gal[True-bad]
    hdrs = [h for h,b in zip(hdrs,bad) if not b]
    log.info("Flagged out %i bad values (%0.1f%%)." % (bad.sum(),bad.sum()/float(bad.size)))

    diagplot(dsub, tsys, noise, dataset, freq=freq, mask=mask, scans=scans,
             **kwargs)
    for xscan in np.unique(obsids):
        match = obsids == xscan
        diagplot(dsub[match], tsys[match], noise[match],
                 dataset+"_obs%i" % xscan, freq=freq, mask=mask, **kwargs)

    return dsub,gal,hdrs


def hdr_to_freq(h):
    freqarr = ((np.arange(h['NCHAN'])+1-h['RCHAN']) * h['FRES'] +
               h['FOFF'] + h['RESTF'])
    return freqarr

def hdr_to_velo(h):
    veloarr = (np.arange(h['NCHAN'])+1-h['RCHAN']) * h['VRES'] + h['VOFF']
    return veloarr

def add_apex_data(data, hdrs, gal, cubefilename, noisecut=np.inf,
                  retfreq=False, excludefitrange=None, varweight=False,
                  debug=False, kernel_fwhm=10./3600.):


    if debug and log.level > 10:
        log.level = 10

    log.info("Data shape: {}".format(data.shape))
    if data.ndim != 2:
        raise ValueError('Data shape is NOT ok.')
    if data.shape[0] != len(hdrs):
        raise ValueError('Data and headers do not match')
    if data.shape[0] != len(gal):
        raise ValueError('Data and coords od not match')

    def data_iterator(data=data, continuum=False, fsw=False):
        shape0 = data.shape[0]
        for ii in xrange(shape0):
        #for ii in xrange(1000):
            yield data[ii,:]

    # as defined on http://www.apex-telescope.org/heterodyne/shfi/het230/lines/
    linefreq = 218222.192
    def velo_iterator(data=None, linefreq=linefreq, headers=hdrs):
        for h in headers:
            if retfreq:
                freqarr = hdr_to_freq(h)
                #veloarr = ((freqarr-linefreq)/linefreq * constants.c).to(u.km/u.s).value
                # needs to be in hz
                yield freqarr*1e6*u.Hz
            else:
                veloarr = hdr_to_velo(h)
                yield veloarr*u.km/u.s

    def coord_iterator(data=None, coordsys_out='galactic', gal=gal):
        for c in gal:
            yield c.l.deg, c.b.deg

    nhits = cubefilename+"_nhits.fits"

    flatheader = fits.getheader(nhits)
    cubeheader = fits.getheader(cubefilename+".fits")

    makecube.add_data_to_cube(cubefilename+".fits", data=data,
                              flatheader=flatheader,
                              cubeheader=cubeheader, linefreq=218.22219,
                              allow_smooth=True,
                              nhits=nhits,
                              data_iterator=data_iterator,
                              coord_iterator=coord_iterator,
                              velo_iterator=velo_iterator,
                              progressbar=True, coordsys='galactic',
                              velocity_offset=0.0, negative_mean_cut=None,
                              add_with_kernel=True, kernel_fwhm=kernel_fwhm,
                              fsw=False,
                              diagnostic_plot_name=None, chmod=False,
                              default_unit=u.GHz if retfreq else u.km/u.s,
                              smoothto=2,
                              noisecut=noisecut,
                              excludefitrange=None,
                              varweight=varweight,
                              continuum_prefix=None)

def add_pipeline_parameters_to_file(fileprefix, pipeline_type, **kwargs):

    if not os.path.exists(fileprefix+".fits"):
        return False

    f = fits.open(fileprefix+".fits")
    f[0].header['PIPECALL'] = (pipeline_type,'build_cube function called')
    for ii,(k,v) in enumerate(kwargs.iteritems()):
        try:
            kw = 'P{pipetype:_<4s}K{n:02d}'.format(n=ii,
                                                   pipetype=pipeline_type)
            keypair = "{k}:{v}".format(k=k, v=v)
            f[0].header[kw] = keypair
        except Exception as ex:
            log.warning("Header could not be updated with key/value pair"
                        "{k}:{v}. Error: {ex}".format(k=k, v=v, ex=ex))
    f.writeto(fileprefix+".fits", clobber=True)

def add_pipeline_header_data(header):
    header['PIPELINE'] = 'Ginsburg 2014 SHFI OTF Pipeline'
    header['TELESCOP'] = 'APEX'
    header['INSTRUME'] = 'SHFI-1'
    header['PIPEDATE'] = (time.strftime("%y_%m_%d_%H:%M:%S"), 'Date pipeline was run')
    from .version import version,githash
    header['PIPEVERS'] = version
    header['PIPEGIT']  = githash
    import sdpy.version
    header['SDPYVERS'] = (sdpy.version.version, 'sdpy version')
    import astropy.version
    header['ASTROPYV'] = (astropy.version.version,'Astropy version')
    try:
        import pyspeckit
        header['PYSPECKV'] = pyspeckit.__version__
    except ImportError,AttributeError:
        pass
    import FITS_tools.version
    header['FITSTOOV'] = (FITS_tools.version.version,'FITS_tools version')
    import scipy.version
    header['SCIPYVER'] = (scipy.version.version,'scipy version')
    import numpy.version
    header['NUMPYVER'] = (numpy.version.version,'numpy version')
    import spectral_cube.version
    header['SPCUBEVE'] = (spectral_cube.version.version,'spectral_cube version')

def make_blanks(gal, header, cubefilename, clobber=True, pixsize=7.2*u.arcsec):

    lrange = (gal.l.wrap_at(180*u.deg).deg.min()+15/3600.,
              gal.l.wrap_at(180*u.deg).deg.max()+15/3600.)
    brange = gal.b.deg.min()+15/3600.,gal.b.deg.max()+15/3600.
    log.info("Map extent automatically determined: "
             "%0.2f < l < %0.2f,  %0.2f < b < %0.2f" % (lrange[0], lrange[1],
                                                        brange[0], brange[1]))

    naxis1 = (lrange[1]-lrange[0])/(pixsize.to(u.deg).value)
    naxis2 = (brange[1]-brange[0])/(pixsize.to(u.deg).value)
    restfreq = (header['RESTF']*u.MHz)
    bmaj = (1.22*10*u.m / restfreq.to(u.m,u.spectral()))*u.radian

    cubeheader, flatheader = makecube.generate_header(np.mean(lrange),
                                                      np.mean(brange),
                                                      naxis1=naxis1,
                                                      naxis2=naxis2,
                                                      naxis3=4096,
                                                      coordsys='galactic',
                                                      ctype3='VRAD',
                                                      bmaj=bmaj.to(u.deg).value,
                                                      bmin=bmaj.to(u.deg).value,
                                                      pixsize=pixsize.to(u.arcsec).value,
                                                      cunit3='km/s',
                                                      output_flatheader='header.txt',
                                                      output_cubeheader='cubeheader.txt',
                                                      cd3=header['VRES'],
                                                      crval3=-1*header['VRES']*header['RCHAN'],
                                                      crpix3=1, clobber=True,
                                                      bunit="K",
                                                      restfreq=restfreq.to(u.Hz).value,
                                                      radio=True)
    add_pipeline_header_data(cubeheader)
    add_pipeline_header_data(flatheader)

    makecube.make_blank_images(cubefilename, cubeheader=cubeheader,
                               flatheader=flatheader, clobber=clobber)

def make_blanks_freq(gal, header, cubefilename, clobber=True, pixsize=7.2*u.arcsec):
    """ complete freq covg """

    lrange = gal.l.wrap_at(180*u.deg).deg.min()+15/3600.,gal.l.wrap_at(180*u.deg).deg.max()+15/3600.
    brange = gal.b.deg.min()+15/3600.,gal.b.deg.max()+15/3600.
    log.info("Map extent: %0.2f < l < %0.2f,  %0.2f < b < %0.2f" % (lrange[0],
                                                                    lrange[1],
                                                                    brange[0],
                                                                    brange[1]))

    naxis1 = int((lrange[1]-lrange[0])/(pixsize.to(u.deg).value)+10)
    naxis2 = int((brange[1]-brange[0])/(pixsize.to(u.deg).value)+10)
    restfreq = (header['RESTF']*u.MHz)
    bmaj = (1.22*10*u.m / restfreq.to(u.m,u.spectral()))*u.radian
    rchan = header['RCHAN']

    #scalefactor = 1./downsample_factor
    #crpix3 = (rchan-1)*scalefactor+0.5+scalefactor/2.

    cubeheader, flatheader = makecube.generate_header(np.mean(lrange),
                                                      np.mean(brange),
                                                      naxis1=naxis1,
                                                      naxis2=naxis2,
                                                      naxis3=header['NCHAN'],
                                                      coordsys='galactic',
                                                      bmaj=bmaj.to(u.deg).value,
                                                      bmin=bmaj.to(u.deg).value,
                                                      pixsize=pixsize.to(u.arcsec).value,
                                                      cunit3='Hz',
                                                      ctype3='FREQ',
                                                      output_flatheader='header.txt',
                                                      output_cubeheader='cubeheader.txt',
                                                      cd3=header['FRES']*1e6,
                                                      crval3=restfreq.to(u.Hz).value,
                                                      crpix3=rchan,
                                                      clobber=True, bunit="K",
                                                      restfreq=restfreq.to(u.Hz).value,
                                                      radio=True)

    add_pipeline_header_data(cubeheader)
    add_pipeline_header_data(flatheader)

    makecube.make_blank_images(cubefilename, flatheader=flatheader,
                               cubeheader=cubeheader, clobber=clobber)


def make_blanks_merge(cubefilename, lowhigh='low', clobber=True,
                      width=1.0*u.GHz, lowest_freq=None, pixsize=7.2*u.arcsec,
                      restfreq=218222.192*u.MHz):
    # total size is 2.3 x 0.4 degrees
    # 1150x
    # center is 0.55 -0.075
    naxis1 = 1150
    naxis2 = 200
    bmaj = (1.22*10*u.m / restfreq.to(u.m,u.spectral()))**-1*u.radian
    cd3 = ((1*u.km/u.s)/constants.c * 218.2*u.GHz).to(u.Hz).value
    naxis3 = int(np.ceil(((width / (218.2*u.GHz) * constants.c) / (u.km/u.s)).decompose().value))
    if lowest_freq is None:
        lowest_freq = 216.8e9 if lowhigh=='low' else 218e9

    cubeheader, flatheader = makecube.generate_header(0.55, -0.075,
                                                      naxis1=naxis1,
                                                      naxis2=naxis2,
                                                      naxis3=naxis3,
                                                      coordsys='galactic',
                                                      bmaj=bmaj.to(u.deg).value,
                                                      bmin=bmaj.to(u.deg).value,
                                                      pixsize=pixsize.to(u.arcsec).value,
                                                      cunit3='Hz',
                                                      ctype3='FREQ',
                                                      output_flatheader='header.txt',
                                                      output_cubeheader='cubeheader.txt',
                                                      cd3=cd3,
                                                      crval3=lowest_freq,
                                                      crpix3=1, clobber=True,
                                                      bunit="K",
                                                      restfreq=restfreq.to(u.Hz).value,
                                                      radio=True)

    add_pipeline_header_data(cubeheader)
    add_pipeline_header_data(flatheader)

    makecube.make_blank_images(cubefilename, flatheader=flatheader,
                               cubeheader=cubeheader, clobber=clobber)

def data_diagplot(data, dataset, ext='png', newfig=False,
                  max_size=1024, freq=None, scans=None,
                  figure=None, axis=None):
    log.info("Doing diagnostics in "+dataset)
    if figure:
        pass
    elif newfig:
        figure = pl.figure()
    else:
        figure = pl.figure(1)
        figure.clf()
    if (np.isnan(data)).all():
        log.exception("ALL data is NaN in {0}".format(dataset))
        import ipdb; ipdb.set_trace()

    if np.any([d > max_size for d in data.shape]):
        # downsample to *not less than* max_size
        factors = [max([1,int(np.floor(d / max_size))]) for d in data.shape]
        data = image_tools.downsample(data, min(factors))

    if axis is None:
        axis = figure.gca()

    axis = mpl_plot_templates.imdiagnostics(data, axis=axis,
                                            second_xaxis=freq)

    if freq is not None:
        #axis.set_xticklabels(np.interp(axis.get_xticks(),
        #                               np.arange(freq.size),
        #                               freq))
        axis.figure.axes[5].set_xlabel("Frequency")
    else:
        axis.set_xlabel("Channel #")
    
    axis.set_ylabel("Integration #")

    if scans is not None and len(scans) < 50:
        xlim = axis.get_xlim()
        ylim = axis.get_ylim()
        axis.hlines(scans, xlim[0], xlim[1], color='k', linestyle='--',
                    alpha=0.5)
        axis.set_xlim(*xlim)

    figfilename = os.path.join(diagplotdir, dataset+"_diagnostics."+ext)
    checkdir_makedir(figfilename)
    try:
        pl.savefig(figfilename,bbox_inches='tight')
    except Exception as ex:
        log.exception(ex)
        print ex
    return axis

def diagplot(data, tsys, noise, dataset, freq=None, mask=None, ext='png',
             newfig=False, **kwargs):
    """
    Generate a set of diagnostic plots

    Parameters
    ----------
    data : `numpy.ndarray`
        A 2D data set, with scans along the y-axis and frequency along the
        x-axis
    tsys : `numpy.ndarray`
        A 1D data set giving TSYS at each time
    noise : `numpy.ndarray`
        The measured noise in each scan
    freq : `numpy.ndarray` or None
        The frequencies to plot along the X-axis
    mask : `numpy.ndarray`
        A boolean mask array with True = good values to be plotted
    ext : str
        The image extension to use when saving
    """

    if newfig:
        pl.figure()
    else:
        pl.figure(2)
        pl.clf()
    pl.subplot(2,1,1)
    pl.plot(tsys,np.arange(tsys.size),alpha=0.5)
    pl.xlabel("TSYS")
    pl.ylabel("Integration")
    pl.subplot(2,1,2)
    pl.plot(tsys, noise, '.',alpha=0.5)
    pl.xlabel("TSYS")
    pl.ylabel("Noise")
    figfilename = os.path.join(diagplotdir, dataset+"_tsys."+ext)
    checkdir_makedir(figfilename)
    pl.savefig(figfilename,bbox_inches='tight')

    if newfig:
        pl.figure()
    else:
        pl.figure(3)
        pl.clf()
    if freq is None:
        freq = np.arange(data.shape[1])
    pl.plot(freq, data.mean(axis=0))
    if mask is not None:
        # Avoid the incorrect appearance of interpolation by masking out
        # intermediate values
        d_to_plot = data.mean(axis=0)
        d_to_plot[mask] = np.nan
        pl.plot(freq, d_to_plot)
    pl.xlabel("Frequency")
    pl.ylabel("Mean Counts")
    figfilename = os.path.join(diagplotdir, dataset+"_masked."+ext)
    checkdir_makedir(figfilename)
    pl.savefig(figfilename,bbox_inches='tight')

    data_diagplot(data, dataset, ext=ext, newfig=newfig, freq=freq, **kwargs)

def build_cube_generic(window, freq=True, mergefile=None, datapath='./',
                       outpath='./', datasets=[], scanblsub=False,
                       shapeselect=None,
                       sourcename=None,
                       tsysrange=[100,250],
                       excludefitrange=None,
                       downsample_factor=None,
                       pixsize=7.2*u.arcsec,
                       kernel_fwhm=10/3600.,
                       pca_clean=False,
                       timewise_pca=True,
                       verbose=False, debug=False, **kwargs):
    """
    TODO: comment!

    kwargs are passed to process_data

    Parameters
    ----------
    window : 'low' or 'high'
        Which of the two APEX SHFI windows to use
    freq : bool
        If True, the cube will be in frequency units and will fully cover the
        observed spectral range.  If False, the cube will be in velocity units
        centered on the observed rest frequency.  This is ignored if mergefile
        is set
    """
    if window not in ('low','high'):
        raise ValueError()
    if mergefile:
        cubefilename=os.path.join(outpath,"{0}_{1}".format(mergefile, window))
    else:
        # assume that we want a cube for EACH data set
        cubefilename = None

    #rcr = [-1000,0] if window == 'low' else [0,5000]
    #xtel = 'AP-H201-F101' if window == 'high' else 'AP-H201-F102'
    xtel = 'AP-H201-X202' if window=='low' else 'AP-H201-X201'

    all_data,all_hdrs,all_gal = {},{},{}
    for dataset in datasets:

        apex_filename = os.path.join(datapath,dataset+".apex")

        spectra,headers,indices = load_apex_cube(apex_filename)
        data,hdrs,gal = select_apex_data(spectra, headers, indices,
                                         sourcename=sourcename,
                                         shapeselect=shapeselect, xtel=xtel,
                                         rchanrange=None,
                                         galactic_coordinate_range=None,
                                         tsysrange=tsysrange)
        log.info("Selected %i spectra from %s" % (len(hdrs), dataset))

        all_data[dataset] = data
        all_hdrs[dataset] = hdrs
        all_gal[dataset] = gal

    all_gal_vect = coordinates.SkyCoord(np.hstack([all_gal[g].l.to(u.radian).value
                                                   for g in all_gal]) * u.radian,
                                        np.hstack([all_gal[g].b.to(u.radian).value
                                                   for g in all_gal]) * u.radian,
                                        frame='galactic')
    all_gal_vect.l.wrap_angle = 180*u.deg

    log.info("Data has been collected and flagged, now adding to cube.")

    headerpars = dict(kernel_fwhm=kernel_fwhm, pca_clean=pca_clean,
                      timewise_pca=timewise_pca,
                      scanblsub=scanblsub)
    if 'pcakwargs' in kwargs:
        headerpars.update(kwargs['pcakwargs'])
    add_pipeline_parameters_to_file(cubefilename, 'generic', **headerpars)
                                    

    for dataset in all_data:

        if not mergefile:
            cubefilename = os.path.join(outpath, "{0}_{1}_cube".format(dataset,window))
            if freq:
                make_blanks_freq(all_gal_vect, hdrs[0], cubefilename,
                                 clobber=True, pixsize=pixsize)
            else:
                make_blanks(all_gal_vect, hdrs[0], cubefilename, clobber=True,
                            pixsize=pixsize)

        data = all_data[dataset]
        hdrs = all_hdrs[dataset]
        gal  = all_gal[dataset]

        data, gal, hdrs = process_data(data, gal, hdrs, dataset+"_"+xtel,
                                       scanblsub=scanblsub, verbose=verbose,
                                       timewise_pca=timewise_pca,
                                       pca_clean=pca_clean, **kwargs)

        add_apex_data(data, hdrs, gal, cubefilename,
                      excludefitrange=excludefitrange,
                      retfreq=freq,
                      varweight=True,
                      kernel_fwhm=kernel_fwhm,
                      debug=debug)

    log.info("Completed cubemaking.  Now doing 'continuum subtraction'"
             "or cube-based spectral baselining")
    cube = fits.open(cubefilename+'.fits', memmap=False)
    cont = fits.getdata(cubefilename+'_continuum.fits')
    data = cube[0].data
    cube[0].data = data - cont
    cube.writeto(cubefilename+'_sub.fits', clobber=True)


    # Downsample by some factor?
    if downsample_factor:
        log.info("Downsampling "+cubefilename)
        avg = FITS_tools.downsample.downsample_axis(cube[0].data, downsample_factor, 0)
        cube[0].data = avg
        cube[0].header['CDELT3'] *= downsample_factor
        scalefactor = 1./downsample_factor
        crpix3 = (cube[0].header['CRPIX3']-1)*scalefactor+0.5+scalefactor/2.
        cube[0].header['CRPIX3'] = crpix3
        cube.writeto(cubefilename+'_downsampled.fits', clobber=True)

    log.info("Done with "+cubefilename)


def build_cube_ao(window, freq=False, mergefile=None,
                  mergefilename=None,
                  datapath=aorawpath,
                  outpath=aopath,
                  datasets=datasets_ao,
                  kernel_fwhm=10/3600.,
                  scanblsub=False,
                  verbose=False,
                  debug=False,
                  pca_clean=True,
                  timewise_pca=True,
                  extra_suffix="",
                  **kwargs):
    """
    TODO: comment!

    kwargs are passed to process_data
    """
    if window not in ('low','high'):
        raise ValueError()
    if mergefile:
        if mergefilename is not None:
            cubefilename = mergefilename
        else:
            cubefilename=os.path.join(outpath,'APEX_H2CO_merge_%s' % window)
    elif freq:
        cubefilename=os.path.join(outpath,'APEX_H2CO_Ao_Freq_%s' % window)
    else:
        cubefilename=os.path.join(outpath,'APEX_H2CO_Ao_%s' % window)

    if extra_suffix:
        cubefilename = cubefilename + extra_suffix

    #rcr = [-1000,0] if window == 'low' else [0,5000]
    xtel = 'AP-H201-F101' if window == 'high' else 'AP-H201-F102'

    all_data,all_hdrs,all_gal = {},{},{}
    for dataset in datasets:

        apex_filename = os.path.join(datapath,dataset+"_merge.apex")

        spectra,headers,indices = load_apex_cube(apex_filename)
        data,hdrs,gal = select_apex_data(spectra, headers, indices,
                                         sourcename='SGRA', shapeselect=4096,
                                         xtel=xtel,
                                         rchanrange=None,
                                         #rchanrange=rcr,
                                         tsysrange=[100,250])
        log.info("Selected %i spectra from %s" % (len(hdrs), dataset))

        #This flagging is more appropriately done in the process_data step
        # # noise_cut = 4 determined by looking at a plot of noise vs time; 0.7%
        # # of data is above 4
        # # Extreme noise appears independent of TSYS!
        # # 4% of data >0.75, but it's pretty bad
        # noise = np.std(data,axis=1)
        # freq_step = np.array([h['FRES'] for h in hdrs])
        # exptime = np.array([h['EXPOSURE'] for h in hdrs])
        # tsys = np.array([h['TSYS'] for h in hdrs])
        # theoretical_rms = 2.0**0.5*tsys/(np.abs(freq_step*1.0e6)*exptime)**0.5
        # bad = noise > theoretical_rms
        # data = data[True-bad]
        # gal = gal[True-bad]
        # hdrs = [h for h,b in zip(hdrs,bad) if not b]
        # print "Flagged out %i bad values (%0.1f%%)." % (bad.sum(),bad.sum()/float(bad.size))

        all_data[dataset] = data
        all_hdrs[dataset] = hdrs
        all_gal[dataset] = gal

    all_gal_vect = coordinates.SkyCoord(np.hstack([all_gal[g].l.to(u.radian).value
                                                   for g in all_gal]) * u.radian,
                                        np.hstack([all_gal[g].b.to(u.radian).value
                                                   for g in all_gal]) * u.radian,
                                        frame='galactic')
    all_gal_vect.l.wrap_angle = 180*u.deg

    if not mergefile:
        if freq:
            make_blanks_freq(all_gal_vect, hdrs[0], cubefilename, clobber=True)
        else:
            make_blanks(all_gal_vect, hdrs[0], cubefilename, clobber=True)

    if freq:
        excludefitrange=None
    else:
        excludefitrange = [700,1300] # FIX THIS when velos are fixed

    headerpars = dict(kernel_fwhm=kernel_fwhm, pca_clean=pca_clean,
                      timewise_pca=timewise_pca,
                      scanblsub=scanblsub)
    if 'pcakwargs' in kwargs:
        headerpars.update(kwargs['pcakwargs'])
    add_pipeline_parameters_to_file(cubefilename, 'ao', **headerpars)

    log.info("Data has been collected and flagged, now adding to cube.")

    for dataset in all_data:
        data = all_data[dataset]
        hdrs = all_hdrs[dataset]
        gal  = all_gal[dataset]

        data, gal, hdrs = process_data(data, gal, hdrs, dataset+"_"+xtel,
                                       scanblsub=scanblsub, verbose=verbose,
                                       pca_clean=pca_clean,
                                       timewise_pca=timewise_pca,
                                       **kwargs)

        add_apex_data(data, hdrs, gal, cubefilename,
                      excludefitrange=excludefitrange,
                      kernel_fwhm=kernel_fwhm,
                      retfreq=freq,
                      varweight=True,
                      debug=debug)


    cube = fits.open(cubefilename+'.fits', memmap=False)
    cont = fits.getdata(cubefilename+'_continuum.fits')
    data = cube[0].data
    cube[0].data = data - cont
    cube.writeto(cubefilename+'_sub.fits', clobber=True)

    if not mergefile:
        # Downsample by averaging over a factor of 8
        downsample_factor = 4 if freq else 8
        avg = np.mean([cube[0].data[ii::downsample_factor,:,:] for ii in
                       xrange(downsample_factor)], axis=0)
        cube[0].data = avg
        cube[0].header['CDELT3'] *= float(downsample_factor)
        scalefactor = 1./downsample_factor
        crpix3 = (cube[0].header['CRPIX3']-1)*scalefactor+0.5+scalefactor/2.
        cube[0].header['CRPIX3'] = crpix3
        # from FITS_tools/hcongrid    h['CRPIX2'] = (h['CRPIX2']-1)*scalefactor + scalefactor/2. + 0.5
        cube.writeto(cubefilename+'_downsampled.fits', clobber=True)

def build_cube_2013(mergefile=None,
                    lowhigh='low',
                    downsample_factor=8,
                    datapath=june2013datapath,
                    outpath=june2013path,
                    datasets=datasets_2013,
                    kernel_fwhm=10/3600.,
                    scanblsub=False,
                    timewise_pca=False, # 2013 data can't handle cleaning.
                    pca_clean=False, # 2013 data can't handle cleaning.  =(
                    extra_suffix="",
                    verbose=True, **kwargs):
    if mergefile:
        cubefilename=os.path.join(outpath,mergefile)
    else:
        cubefilename=os.path.join(outpath,
                                  'APEX_H2CO_2013_%s' % lowhigh)
    if extra_suffix:
        cubefilename = cubefilename + extra_suffix

    xtel = 'AP-H201-X202' if lowhigh=='low' else 'AP-H201-X201'

    if not mergefile:
        # Need two loops.  First one is just to determine map extent.
        all_gal = {}
        for dataset in datasets:

            apex_filename=datapath+dataset+".apex"

            spectra,headers,indices = load_apex_cube(apex_filename,
                                                     skip_data=True,
                                                     downsample_factor=downsample_factor)
            data, hdrs, gal = select_apex_data(spectra, headers, indices,
                                               sourcename='SGRA',
                                               shapeselect=32768/downsample_factor,
                                               tsysrange=[100,325],
                                               xtel=xtel,
                                               rchanrange=None,
                                               skip_data=True)
            all_gal[dataset] = gal

        all_gal_vect = coordinates.SkyCoord(np.hstack([all_gal[g].l.to(u.radian).value
                                                       for g in all_gal]) * u.radian,
                                            np.hstack([all_gal[g].b.to(u.radian).value
                                                       for g in all_gal]) * u.radian,
                                            frame='galactic')
        all_gal_vect.l.wrap_angle = 180*u.deg

        make_blanks_freq(all_gal_vect, hdrs[0], cubefilename, clobber=True)

    headerpars = dict(kernel_fwhm=kernel_fwhm, pca_clean=pca_clean,
                      timewise_pca=timewise_pca,
                      scanblsub=scanblsub)
    if 'pcakwargs' in kwargs:
        headerpars.update(kwargs['pcakwargs'])
    add_pipeline_parameters_to_file(cubefilename, '2013', **headerpars)

    # need two loops to avoid loading too much stuff into memory
    for dataset in datasets:

        log.info("Adding data set {0} to cube file {1}".format(dataset, cubefilename))

        apex_filename=datapath+dataset+".apex"

        spectra,headers,indices = load_apex_cube(apex_filename, skip_data=False,
                                                 downsample_factor=downsample_factor,
                                                 )

        if dataset == 'M-091.F-0019-2013-2013-06-13':
            tsysrange=[100,260]
        else:
            tsysrange=[100,325]

        data, hdrs, gal = select_apex_data(spectra, headers, indices,
                                           sourcename='SGRA',
                                           # NOT ignored, even though it's not used above...
                                           # this is only OK because the bad shapes are from
                                           # Saturn
                                           shapeselect=32768/downsample_factor,
                                           tsysrange=tsysrange,
                                           xtel=xtel,
                                           rchanrange=None,
                                           skip_data=False)
        
        data, gal, hdrs = process_data(data, gal, hdrs, dataset+"_"+xtel,
                                       scanblsub=scanblsub, verbose=verbose,
                                       timewise_pca=timewise_pca,
                                       pca_clean=pca_clean,
                                       **kwargs)

        add_apex_data(data, hdrs, gal, cubefilename, retfreq=True,
                      kernel_fwhm=kernel_fwhm,
                      varweight=True,)
        # FORCE cleanup
        del data,hdrs,gal

    cube = fits.open(cubefilename+'.fits', memmap=False)
    cont = fits.getdata(cubefilename+'_continuum.fits')
    data = cube[0].data
    cube[0].data = data - cont
    cube.writeto(cubefilename+'_sub.fits', clobber=True)

    # Downsample by averaging over a factor of 8
    # (this is extra downsampling)
    avg = np.mean([cube[0].data[ii::2,:,:] for ii in xrange(2)], axis=0)
    cube[0].data = avg
    cube[0].header['CDELT3'] *= 2
    scalefactor = 1./2.
    crpix3 = (cube[0].header['CRPIX3']-1)*scalefactor+0.5+scalefactor/2.
    cube[0].header['CRPIX3'] = crpix3
    cube.writeto(cubefilename+'_downsampled.fits', clobber=True)

def build_cube_2014(sourcename,
                    mergefile=None,
                    lowhigh='low',
                    downsample_factor=8,
                    datapath=april2014path,
                    kernel_fwhm=10/3600.,
                    outpath=april2014path,
                    datasets=None,
                    scanblsub=False,
                    verbose=True,
                    pca_clean=False,
                    timewise_pca=False,
                    extra_suffix='',
                    **kwargs
                    ):
    """
    Wrapper.  Because each field has its own name in 2014, this will need to be
    modified for the mergefile to accept wildcards or something for sourcename
    selection
    """
    if mergefile:
        cubefilename=os.path.join(outpath,mergefile)
    elif isinstance(sourcename, str):
        cubefilename=os.path.join(outpath, 
                                  'APEX_H2CO_2014_%s_%s' % (sourcename, lowhigh))
    else:
        raise ValueError("Use a mergefile")
    if extra_suffix:
        cubefilename = cubefilename + extra_suffix

    log.info("Building cubes for "+cubefilename)

    xtel = 'AP-H201-X202' if lowhigh=='low' else 'AP-H201-X201'

    t0 = time.time()

    if not mergefile:
        # Need two loops.  First one is just to determine map extent.
        all_gal = {}
        for dataset in datasets:

            apex_filename=datapath+dataset+".apex"

            log.info("".join(("Pre-Loading data for dataset ",dataset," to filename ",apex_filename,"  t=",str(time.time()-t0))))

            spectra,headers,indices = load_apex_cube(apex_filename,
                                                     skip_data=True,
                                                     downsample_factor=downsample_factor)
            data, hdrs, gal = select_apex_data(spectra, headers, indices,
                                               sourcename=sourcename,
                                               shapeselect=32768/downsample_factor,
                                               tsysrange=[100,325],
                                               xtel=xtel,
                                               rchanrange=None,
                                               skip_data=True)
            all_gal[dataset] = gal

        all_gal_vect = coordinates.SkyCoord(np.hstack([all_gal[g].l.to(u.radian).value
                                                       for g in all_gal]) * u.radian,
                                            np.hstack([all_gal[g].b.to(u.radian).value
                                                       for g in all_gal]) * u.radian,
                                            frame='galactic')
        all_gal_vect.l.wrap_angle = 180*u.deg

        log.info("Making blanks for "+cubefilename)
        make_blanks_freq(all_gal_vect, hdrs[0], cubefilename, clobber=True)

    headerpars = dict(kernel_fwhm=kernel_fwhm, pca_clean=pca_clean,
                      timewise_pca=timewise_pca,
                      scanblsub=scanblsub)
    if 'pcakwargs' in kwargs:
        headerpars.update(kwargs['pcakwargs'])
    add_pipeline_parameters_to_file(cubefilename, '2014', **headerpars)

    # need two loops to avoid loading too much stuff into memory
    for dataset in datasets:

        apex_filename=datapath+dataset+".apex"
        
        log.info("".join(("Loading data for dataset ",dataset," in filename ",apex_filename,"  t=",str(time.time()-t0))))

        spectra,headers,indices = load_apex_cube(apex_filename, skip_data=False,
                                                 downsample_factor=downsample_factor,
                                                 )

        #if dataset == 'M-091.F-0019-2013-2013-06-13':
        #    tsysrange=[100,260]
        #else:
        #    tsysrange=[100,325]
        tsysrange=[100,325]

        log.info("".join(("Selecting data for dataset ",dataset," in filename ",apex_filename,"  t=",str(time.time()-t0))))

        data, hdrs, gal = select_apex_data(spectra, headers, indices,
                                           sourcename=sourcename,
                                           # NOT ignored, even though it's not used above...
                                           # this is only OK because the bad shapes are from
                                           # Saturn
                                           #shapeselect=4096,
                                           shapeselect=32768/downsample_factor,
                                           tsysrange=tsysrange,
                                           xtel=xtel,
                                           rchanrange=None,
                                           skip_data=False)

        log.info("".join(("Processing data for dataset ",dataset," in filename ",apex_filename,"  t=",str(time.time()-t0))))

        data, gal, hdrs = process_data(data, gal, hdrs, os.path.join(outpath,
                                                                     dataset)+"_"+xtel,
                                       scanblsub=scanblsub, verbose=verbose,
                                       timewise_pca=timewise_pca,
                                       pca_clean=pca_clean,
                                       **kwargs)

        log.info("".join(("Adding data for dataset ",dataset," to filename ",cubefilename,"  t=",str(time.time()-t0))))

        add_apex_data(data, hdrs, gal, cubefilename, retfreq=True,
                      kernel_fwhm=kernel_fwhm,
                      varweight=True,
                      # downsample factor for freqarr
                      )
        # FORCE cleanup
        log.info("".join(("Clearing data for dataset ",dataset," to filename ",cubefilename,"  t=",str(time.time()-t0))))
        del data,hdrs,gal

    log.info("".join(("Continuum subtracting ",cubefilename)))

    cube = fits.open(cubefilename+'.fits', memmap=False)
    cont = fits.getdata(cubefilename+'_continuum.fits')
    data = cube[0].data
    cube[0].data = data - cont
    cube.writeto(cubefilename+'_sub.fits', clobber=True)

    log.info("Downsampling "+cubefilename)

    # Downsample by averaging over a factor of 8
    avg = FITS_tools.downsample.downsample_axis(cube[0].data, 2, 0)
    cube[0].data = avg
    cube[0].header['CDELT3'] *= 2
    scalefactor = 1./2.
    crpix3 = (cube[0].header['CRPIX3']-1)*scalefactor+0.5+scalefactor/2.
    cube[0].header['CRPIX3'] = crpix3
    cube.writeto(cubefilename+'_downsampled.fits', clobber=True)

    log.info("Done with "+cubefilename)




def make_high_mergecube(pca_clean={'2014':True,
                                   '2013':False,
                                   'ao':True},
                        scanblsub={'2014':False, '2013':False, 'ao':False},
                        timewise_pca={'2014': True, '2013':False, 'ao':True},
                        mergefile2=None):

    if mergefile2 is None:
        raise ValueError("Must specify a merge filename")
    #if pca_clean:
    #    if timewise_pca:
    #        mergefile2 = 'APEX_H2CO_merge_high_timepca'
    #    else:
    #        mergefile2 = 'APEX_H2CO_merge_high'
    #else:
    #    mergefile2 = 'APEX_H2CO_merge_high_nopca'
    make_blanks_merge(os.path.join(mergepath,mergefile2), lowhigh='high',
                      lowest_freq=218e9, width=1.0*u.GHz)

    mapnames = ['MAP_{0:03d}'.format(ii) for ii in range(1,130)]
    log.info("Building cubes: "+str(mapnames))
    # Frequency: (216.9, 219.4)
    build_cube_2014(mapnames,
                    mergefile=mergefile2,
                    outpath=mergepath,
                    datapath=april2014path,
                    lowhigh='low',
                    pca_clean=pca_clean['2014'],
                    timewise_pca=timewise_pca['2014'],
                    scanblsub=scanblsub['2014'],
                    datasets=datasets_2014)

    log.info("Building Ao cubes")
    # ('ao', 'high'): (218.0, 219.0),
    build_cube_ao(window='high', mergefile=True, freq=True, outpath=mergepath,
                  pca_clean=pca_clean['ao'], timewise_pca=timewise_pca['ao'],
                  mergefilename=os.path.join(mergepath, mergefile2),
                  scanblsub=scanblsub['ao'],
                  datapath=aorawpath)

    log.info("Building 2013 cubes")
    # (2013, 'high'): (217.5, 220.0)
    build_cube_2013(mergefile=mergefile2,
                    outpath=mergepath,
                    datapath=june2013datapath,
                    lowhigh='high',
                    timewise_pca=timewise_pca['2013'],
                    pca_clean=pca_clean['2013'],
                    scanblsub=scanblsub['2013'])



def make_low_mergecube(datasets_2014=datasets_2014, pca_clean=True,
                       scanblsub=False, timewise_pca=True, mergefile1=None):
    mergefile1 = 'APEX_H2CO_merge_low'
    make_blanks_merge(os.path.join(mergepath,mergefile1), lowhigh='low')
    build_cube_ao(window='low',
                  mergefile=True,
                  freq=True,
                  outpath=mergepath,
                  datapath=june2013datapath,
                  timewise_pca=timewise_pca,
                  pca_clean=pca_clean,
                  scanblsub=scanblsub,
                  mergefilename=os.path.join(mergepath, mergefile1),
                 )
    build_cube_2013(mergefile=mergefile1,
                    outpath=mergepath,
                    datapath=june2013datapath,
                    lowhigh='low',
                    scanblsub=scanblsub)

    log.info("Starting merge")
    for lowhigh in ('low',):
        mapnames = ['MAP_{0:03d}'.format(ii) for ii in range(1,130)]
        log.info("Building cubes: "+str(mapnames)+" "+lowhigh)
        build_cube_2014(mapnames,
                        mergefile=mergefile1,
                        outpath=mergepath,
                        lowhigh=lowhigh,
                        datasets=datasets_2014)


def integrate_slices_high(prefix='merged_datasets/APEX_H2CO_merge_high_sub'):
    ffile = fits.open(prefix+'.fits')

    integ1,hdr = cubes.integ(ffile, [235,344], average=np.nansum) # first H2CO line: blue
    hdu1 = fits.PrimaryHDU(data=integ1, header=hdr)
    hdu1.writeto(prefix+"_H2CO_303-202_blue.fits", clobber=True)
    integ2,hdr = cubes.integ(ffile, [161,235], average=np.nansum) # first H2CO line: red
    hdu2 = fits.PrimaryHDU(data=integ2, header=hdr)
    hdu2.writeto(prefix+"_H2CO_303-202_red.fits", clobber=True)


    integ4,hdr = cubes.integ(ffile, [161,344], average=np.nansum) # first H2CO line: red
    hdu4 = fits.PrimaryHDU(data=integ4, header=hdr)
    hdu4.writeto(prefix+"_H2CO_303-202.fits", clobber=True)


    integ3,hdr = cubes.integ(ffile, [513,615], average=np.nansum) # second H2CO line: blue
    hdu3 = fits.PrimaryHDU(data=integ3, header=hdr)
    hdu3.writeto(prefix+"_H2CO_322-221_blue.fits", clobber=True)

def integrate_slices_low(prefix='merged_datasets/APEX_H2CO_merge_low_sub'):
    ffile = fits.open(prefix+'.fits')

    integ1,hdr = cubes.integ(ffile, [335,446], average=np.nansum)
    hdu1 = fits.PrimaryHDU(data=integ1, header=hdr)
    hdu1.writeto(prefix+"_SiO5-4.fits", clobber=True)

def integrate_mask(prefix, mask=h2copath+'APEX_H2CO_303_202_mask.fits'):
    """
    Integrate a cube with name specified by 'prefix' using a specific mask
    """
    if isinstance(mask,str):
        mask = fits.getdata(mask).astype('bool')
    ffile = fits.open(prefix+'.fits')
    cd = ffile[0].header['CDELT3']
    ffile[0].data *= mask * cd
    ffile[0].data[~mask.astype('bool')] = np.nan

    integ1,hdr = cubes.integ(ffile, [0,ffile[0].shape[0]], average=np.nansum)
    hdr['BUNIT'] = ('K km/s',"Integrated over masked region")
    hdu1 = fits.PrimaryHDU(data=integ1, header=hdr)
    hdu1.writeto(prefix+"_mask_integ.fits", clobber=True)

def integrate_h2co_by_freq(filename):
    import spectral_cube
    cube = spectral_cube.SpectralCube.read(filename)

    #if 'high' in filename:
    #    cocube = cube
    #else:
    #    cocube = spectral_cube.SpectralCube.read(filename.replace('low','high'))

    #mcube = cocube.with_spectral_unit(u.km/u.s,
    #                                rest_value=bright_lines['13CO']*u.GHz,
    #                                velocity_convention='radio')
    #coscube = mcube.spectral_slab(-100*u.km/u.s, 150*u.km/u.s)
    #mask = coscube > 1

    for line in bright_lines:
        scube = cube.with_spectral_unit(u.km/u.s,
                                        rest_value=bright_lines[line]*u.GHz,
                                        velocity_convention='radio')
        subcube1 = scube.spectral_slab(-100*u.km/u.s, 150*u.km/u.s)
        ncube = scube.spectral_slab(-150*u.km/u.s, -100*u.km/u.s)
        noise = ncube.apply_numpy_function(np.std, axis=0)
        #mask._wcs = subcube1.wcs
        subcube = subcube1.with_mask(subcube1>noise)#.with_mask(mask)
        if subcube.shape[0] == 1:
            # implies out of range
            continue
        mom0 = subcube.moment0()
        mom1 = subcube.moment1()
        mom2 = subcube.moment2()
        fn = os.path.split(filename)[1]
        outfn = 'projections/'+fn.replace(".fits","_{line}_{mom}.fits")
        mom0.hdu.writeto(outfn.format(line=line, mom='mom0'),clobber=True)
        mom1.hdu.writeto(outfn.format(line=line, mom='mom1'),clobber=True)
        mom2.hdu.writeto(outfn.format(line=line, mom='mom2'),clobber=True)

def compute_noise_high(prefix=mergepath+'APEX_H2CO_merge_high_sub',pixrange=[700,900]):
    ffile = fits.open(prefix+'.fits')

    integ1,hdr = cubes.integ(ffile, pixrange, average=np.nanstd)
    hdu1 = fits.PrimaryHDU(data=integ1, header=hdr)
    hdu1.writeto(prefix+"_noise.fits", clobber=True)

def compute_noise_low(prefix=mergepath+'APEX_H2CO_merge_low_sub',pixrange=[512,675]):
    ffile = fits.open(prefix+'.fits')

    integ1,hdr = cubes.integ(ffile, pixrange, average=np.nanstd)
    hdu1 = fits.PrimaryHDU(data=integ1, header=hdr)
    hdu1.writeto(prefix+"_noise.fits", clobber=True)

def compute_noise_extras(prefix=june2013path+'APEX_H2CO_2013_%s_sub',
                         lowhigh='high',
                         pixrange=[0,4096]):
    ffile = fits.open((prefix % lowhigh)+'.fits')

    integ1,hdr = cubes.integ(ffile, pixrange, average=np.nanstd)
    hdu1 = fits.PrimaryHDU(data=integ1, header=hdr)
    hdu1.writeto(prefix+"_noise.fits", clobber=True)

def signal_to_noise_mask_cube(prefix, noise=None, kernelsize=[2,2,2], grow=1,
                              sigmacut=3):
    """
    It's not clear that this has anything to do with moments...
    (this was called moment_mask_cube in the past, which I think is just wrong)
    """
    ffile = fits.open(prefix+'.fits')
    cube = ffile[0].data
    if noise is None:
        noise = fits.getdata(prefix+'_noise.fits')

    t0 = time.time()
    log.info("Initiating cube smooth of {0}.".format(prefix))
    smcube = cube_regrid.gsmooth_cube(cube, kernelsize, use_fft=False,
                                      kernelsize_mult=3)
    log.info("Completed cube smooth in %i seconds" % (time.time()-t0))
    mask = smcube > noise*sigmacut

    mask_grow = scipy.ndimage.morphology.binary_dilation(mask, iterations=grow)

    ffile[0].data[True-mask_grow] = np.nan
    ffile[0].writeto(prefix+"_snmasked.fits", clobber=True)

    ffile[0].data = mask_grow.astype('int')
    ffile[0].writeto(prefix+"_mask.fits", clobber=True)

def do_sncube_masking_hi(prefix=h2copath+'APEX_H2CO_303_202'):
    # 0-25 not checked! arbitrary choice.
    compute_noise_high(prefix, pixrange=[0,25])
    signal_to_noise_mask_cube(prefix)
    integrate_slices_high(prefix+'_snmasked')

def extract_subcube(cubefilename, outfilename, linefreq=218.22219*u.GHz,
                    debug=False, smooth=False, vsmooth=False, naxis3=400,
                    vmin=-150*u.km/u.s, vmax=250*u.km/u.s):
    t0 = time.time()
    log.info(("Extracting subcube at {0} from {1}"
              " with smooth={2} and vsmooth={3}").format(linefreq,
                                                         cubefilename, smooth,
                                                         vsmooth))

    cube = spectral_cube.SpectralCube.read(cubefilename)
    vcube = cube.with_spectral_unit(u.km/u.s, rest_value=linefreq,
                                    velocity_convention='radio')
    svcube = vcube.spectral_slab(vmin, vmax)
    crval3 = vmin.to(u.km/u.s).value

    outheader = svcube.header
    outheader['CRPIX3'] = 1
    outheader['CRVAL3'] = crval3
    outheader['CUNIT3'] = 'km/s'
    outheader['CDELT3'] = 1.0
    outheader['NAXIS3'] = naxis3
    outheader['NAXIS2'] = svcube.shape[1]
    outheader['NAXIS1'] = svcube.shape[2]

    if smooth:
        #cubesm = gsmooth_cube(ffile[0].data, [3,2,2], use_fft=True,
        #                      psf_pad=False, fft_pad=False)
        # smoothed with 2 pixels -> sigma=10", fwhm=23"
        # this is an "optimal smooth", boosting s/n and smoothing to 36"
        # resolution.
        kw = 2 if not vsmooth else 4
        cubesm = cube_regrid.spatial_smooth_cube(svcube.filled_data[:], kw,
                                                 use_fft=False,
                                                 numcores=4)
        cubesm = cube_regrid.spectral_smooth_cube(cubesm, 3/2.35,
                                                  use_fft=False,
                                                  numcores=4)
        svcube._data = cubesm

        outheader['CDELT3'] = outheader['CDELT3'] * kw
        outheader['NAXIS3'] = outheader['NAXIS3'] / kw
    
    # Now that we've written this out, we use interpolation to force the cube
    # onto a grid that starts at *exactly* -150 km/s
    newhdu = cube_regrid.regrid_cube_hdu(svcube.hdu, outheader, order=1, prefilter=False)
    newhdu.writeto(outfilename, output_verify='fix', clobber=True)

    log.info("Completed cube extraction to {1} in {0} seconds.".format(time.time()-t0,
                                                                       outfilename))
    
    return newhdu

all_lines = {'H2CO_303_202':218.22219,
             'H2CO_322_221':218.47563,
             'H2CO_321_220':218.76007,
             'SiO_54':217.10498,
             'CH3OH_422_312':218.44005,
             'CH3OH_514_422':216.9456,
             'CH3OH_633_716':216.85786,
             'HCCCH_65': 217.82215,
             'OCS_18_17':218.90336,
             'CH3OCHO_17_16':218.29789,
             'C18O':219.56036,
             '13CO':220.39868,
             #'H2S 2(2,0)-2(1,1)': 216.71044, ??
             }
bright_lines = {k:all_lines[k] for k in
                ['H2CO_303_202', 'H2CO_322_221', 'H2CO_321_220', 'SiO_54',
                 'CH3OH_422_312', 'CH3OH_514_422', 'CH3OH_633_716', 'C18O',
                 '13CO']}
bandwidths = {'H2CO_303_202':25,
              'H2CO_322_221':25,
              'H2CO_321_220':25,
              'SiO_54':25,
              'CH3OH_422_312':25,
              'CH3OH_514_422':25,
              'CH3OH_633_716':25,
              'HCCCH_65': 25,
              'OCS_18_17':25,
              'CH3OCHO_17_16':25,
              'C18O':75,
              '13CO':75,
              #'H2S 2(2,0)-2(1,1)': 216.71044, ??
              }

lines218 = {x:v for x,v in all_lines.iteritems()
            if 'H2CO' in x or 'CH3OH_422_312' in x}

def make_line_mask(freqarr, lines=bright_lines):
    mask = np.ones(freqarr.size, dtype='bool')
    for ln,lf in lines.iteritems():
        bw = bandwidths[ln]
        wh = (lf*1e3-bw < freqarr) & (lf*1e3+bw > freqarr)
        mask[wh] = False
    return mask


def do_extract_subcubes(outdir=mergepath, merge_prefix='APEX_H2CO_merge',
                        frange=None, lines=all_lines):

    for line,freq in lines.iteritems():
        if frange is not None:
            if freq<frange[0] or freq>frange[1]:
                log.info("Skipping line {0}".format(line))
                continue

        log.info("Extracting {0}".format(line))
        if freq < 218:
            cubefilename = os.path.join(mergepath, merge_prefix+'_low_sub.fits')
        else:
            cubefilename = os.path.join(mergepath, merge_prefix+'_high_sub.fits')

        header = fits.getheader(cubefilename)
        ww = wcs.WCS(header)
        wspec = ww.sub([wcs.WCSSUB_SPECTRAL])
        nax = header['NAXIS%i' % (ww.wcs.spec+1)]
        freqarr = wspec.wcs_pix2world(np.arange(nax),0)[0]
        if freq*1e9 > freqarr.min() and freq*1e9 < freqarr.max():
            extract_subcube(cubefilename,
                            os.path.join(outdir, 'APEX_{0}_smooth.fits').format(line),
                            linefreq=freq*u.GHz, smooth=True)
            extract_subcube(cubefilename,
                            os.path.join(outdir, 'APEX_{0}_vsmooth.fits').format(line),
                            linefreq=freq*u.GHz, smooth=True, vsmooth=True)
            extract_subcube(cubefilename,
                            os.path.join(outdir, 'APEX_{0}.fits').format(line),
                            linefreq=freq*u.GHz)
        else:
            log.info("Skipping line {0}".format(line))


def do_everything(pca_clean={'2014':False, '2013':False, 'ao':False},
                  scanblsub={'2014':False, '2013':False, 'ao':False},
                  timewise_pca={'2014':True, '2013':False, 'ao':True},
                  mergefile2='APEX_H2CO_merge_high'):
    make_high_mergecube(mergefile2=mergefile2, pca_clean=pca_clean,
                        scanblsub=scanblsub, timewise_pca=timewise_pca)

    do_postprocessing()


def do_postprocessing():
    #make_low_mergecube() # there's only one really useful overlap region
    os.chdir(mergepath)
    # vsmoothds is made here:
    os.system('./APEX_H2CO_merge_high_starlink_custom.sh')
    os.chdir('../')
    do_extract_subcubes(outdir=mergepath, frange=[218,219],
                        lines=lines218)
    compute_noise_high(mergepath+'APEX_H2CO_merge_high_sub', pixrange=[700,900])
    compute_noise_high(mergepath+'APEX_H2CO_merge_high_smooth',[203,272])
    compute_noise_high(mergepath+'APEX_H2CO_merge_high_vsmoothds',[203,272])
    compute_noise_high(mergepath+'APEX_H2CO_303_202_vsmooth',[75,100])
    #compute_noise_low()
    signal_to_noise_mask_cube(mergepath+'APEX_H2CO_303_202',
                              noise=fits.getdata(mergepath+'APEX_H2CO_merge_high_sub_noise.fits'),
                              sigmacut=2,
                              grow=2)
    signal_to_noise_mask_cube(mergepath+'APEX_H2CO_303_202_smooth',
                              noise=fits.getdata(mergepath+'APEX_H2CO_merge_high_smooth_noise.fits'),
                              sigmacut=2)
    signal_to_noise_mask_cube(mergepath+'APEX_H2CO_303_202_vsmooth',
                              noise=fits.getdata(mergepath+'APEX_H2CO_303_202_vsmooth_noise.fits'),
                              sigmacut=2)
    integrate_mask(mergepath+'APEX_H2CO_303_202',
                   mask=mergepath+'APEX_H2CO_303_202_mask.fits')
    integrate_mask(mergepath+'APEX_H2CO_303_202_smooth',
                   mask=mergepath+'APEX_H2CO_303_202_smooth_mask.fits')
    integrate_mask(mergepath+'APEX_H2CO_303_202_vsmooth',
                   mask=mergepath+'APEX_H2CO_303_202_vsmooth_mask.fits')

    for fn in glob.glob(os.path.join(mergepath,'APEX_H2CO_30*fits')):
        try:
            os.symlink(fn,
                       os.path.join(h2copath,os.path.split(fn)[-1]))
        except OSError:
            log.debug("Skipped file {0} because it exists".format(fn))

    # Create a few integrated H2CO 303 maps
    do_sncube_masking_hi()

    # Use spectral_cube to do a bunch of integrations
    # PATH SENSITIVE
    # integrate_h2co_by_freq(mergepath+mergefile2+".fits")
    # On second thought, let's not go to camelot
    # (this function proved ineffective)

    for line in lines218:
        fn = mergepath+'APEX_{0}.fits'.format(line)
        if os.path.exists(fn):
            integrate_mask(mergepath+'APEX_{0}'.format(line),
                           mask=mergepath+'APEX_H2CO_303_202_mask.fits')
            integrate_mask(mergepath+'APEX_{0}_smooth'.format(line),
                           mask=mergepath+'APEX_H2CO_303_202_smooth_mask.fits')
            integrate_mask(mergepath+'APEX_{0}_vsmooth'.format(line),
                           mask=mergepath+'APEX_H2CO_303_202_vsmooth_mask.fits')
            log.debug("Integrated masked file {0}".format(fn))
        else:
            log.debug("File {0} does not exist".format(fn))

    for line in lines218:
        if os.path.exists(mergepath+'APEX_{0}.fits'.format(line)):
            baseline_cube(mergepath+'APEX_{0}.fits'.format(line),
                          maskfn=mergepath+'APEX_H2CO_303_202_smooth_mask.fits')
            baseline_cube(mergepath+'APEX_{0}_smooth.fits'.format(line),
                          maskfn=mergepath+'APEX_H2CO_303_202_smooth_mask.fits')
            baseline_cube(mergepath+'APEX_{0}_vsmooth.fits'.format(line),
                          maskfn=mergepath+'APEX_H2CO_303_202_vsmooth_mask.fits')

    #compute_noise_high(mergepath+'APEX_H2CO_303_202_bl',[350,400])
    #compute_noise_high(mergepath+'APEX_H2CO_303_202_smooth_bl',[175,200])
    #compute_noise_high(mergepath+'APEX_H2CO_303_202_vsmooth_bl',[80,100])
    signal_to_noise_mask_cube(mergepath+'APEX_H2CO_303_202_bl',
                              noise=fits.getdata(mergepath+'APEX_H2CO_merge_high_sub_noise.fits'),
                              grow=2)
    signal_to_noise_mask_cube(mergepath+'APEX_H2CO_303_202_smooth_bl',
                              noise=fits.getdata(mergepath+'APEX_H2CO_merge_high_smooth_noise.fits'),
                              sigmacut=3)
    signal_to_noise_mask_cube(mergepath+'APEX_H2CO_303_202_vsmooth_bl',
                              noise=fits.getdata(mergepath+'APEX_H2CO_303_202_vsmooth_noise.fits'),
                              sigmacut=3)

    for line in lines218:
        if os.path.exists(mergepath+'APEX_{0}_bl.fits'.format(line)):
            integrate_mask(mergepath+'APEX_{0}_bl'.format(line),
                           mask=mergepath+'APEX_H2CO_303_202_mask.fits')
            integrate_mask(mergepath+'APEX_{0}_smooth_bl'.format(line),
                           mask=mergepath+'APEX_H2CO_303_202_smooth_mask.fits')
            integrate_mask(mergepath+'APEX_{0}_vsmooth_bl'.format(line),
                           mask=mergepath+'APEX_H2CO_303_202_vsmooth_mask.fits')

    do_mask_ch3oh()

    for fn in glob.glob(os.path.join(mergepath,'APEX_H2CO_3*fits')):
        try:
            os.symlink(fn,
                       os.path.join(h2copath,os.path.split(fn)[-1]))
            log.info("Linked file {0} to {1}".format(fn, h2copath))
        except OSError:
            log.debug("Skipped file {0} because it exists".format(fn))

    do_temperature()

def baseline_cube(cubefn, maskfn, order=5):
    from pyspeckit.cubes.cubes import baseline_cube
    f = fits.open(cubefn)
    cube = f[0].data
    mask = fits.getdata(maskfn).astype('bool')
    t0 = time.time()
    log.info("Baselining cube {0} with order {1}...".format(cubefn, order))
    bc = baseline_cube(cube, order, cubemask=mask)
    log.info("Baselining done ({0} seconds)".format(time.time()-t0))
    f[0].data = bc
    f.writeto(cubefn.replace(".fits","_bl.fits"), clobber=True)


def do_everything_2013extrafreqs():
    build_cube_2013(lowhigh='low',
                    scanblsub=False)
    build_cube_2013(lowhigh='high',
                    scanblsub=False)
    #raise NotImplementedError
    #compute_noise_extras(lowhigh='low',pixrange=[0,4096])
    #compute_noise_extras(lowhigh='high',pixrange=[0,4096])


def doratio(h2copath=h2copath, maxratio=1):
    """
    I swapped top and bottom because that's what the models were set up for... 

    maxratio : float
        0.64 is the highest physical value
    
    """
    ##### integrated #####

    for smooth in ('_bl','_smooth_bl','_vsmooth_bl'):
        top = fits.getdata(h2copath+'APEX_H2CO_303_202{0}_mask_integ.fits'.format(smooth))
        bottom = fits.getdata(h2copath+'APEX_H2CO_322_221{0}_CH3OHchomped_mask_integ.fits'.format(smooth))
        
        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>maxratio] = np.nan
        
        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask_integ.fits'.format(smooth))
        
        f[0].data = ratio
        
        f.writeto(h2copath+'H2CO_322221_to_303202{0}_integ.fits'.format(smooth),clobber=True)

        bottom = fits.getdata(h2copath+'APEX_H2CO_321_220{0}_mask_integ.fits'.format(smooth))
        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>maxratio] = np.nan
        
        f[0].data = ratio
        
        f.writeto(h2copath+'H2CO_321220_to_303202{0}_integ.fits'.format(smooth),clobber=True)

    ##### cube #####
    for smooth in ('_bl','_smooth_bl','_vsmooth_bl'):
        top = fits.getdata(h2copath+'APEX_H2CO_303_202{0}.fits'.format(smooth))
        bottom = fits.getdata(h2copath+'APEX_H2CO_322_221{0}_CH3OHchomped_masked.fits'.format(smooth))
        
        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>maxratio] = np.nan
        
        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask.fits'.format(smooth))
        
        # mask now!
        f[0].data = ratio * f[0].data.astype('bool')
        
        f.writeto(h2copath+'H2CO_322221_to_303202_cube{0}.fits'.format(smooth),clobber=True)

        bottom = fits.getdata(h2copath+'APEX_H2CO_321_220{0}.fits'.format(smooth))
        
        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>maxratio] = np.nan
        
        #f = fits.open(h2copath+'APEX_H2CO_303_202{0}.fits'.format(smooth))
        f = fits.open(h2copath+'APEX_H2CO_303_202{0}_mask.fits'.format(smooth))
        
        # mask out...
        f[0].data = ratio * f[0].data.astype('bool')
        
        f.writeto(h2copath+'H2CO_321220_to_303202_cube{}.fits'.format(smooth),clobber=True)



def dopeaksn():

    from FITS_tools import strip_headers

    f = fits.open(h2copath+'APEX_H2CO_303_202.fits')
    header = strip_headers.flatten_header(f[0].header)
    f[0].header=header
    f[0].data = f[0].data.max(axis=0)
    n = fits.getdata(h2copath+'APEX_H2CO_merge_high_sub_noise.fits')
    f[0].data /= n
    f.writeto(h2copath+'APEX_H2CO_303_202_peaksn.fits',clobber=True)

    f = fits.open(h2copath+'APEX_H2CO_303_202_smooth.fits')
    header = strip_headers.flatten_header(f[0].header)
    f[0].header=header
    f[0].data = f[0].data.max(axis=0)
    n = fits.getdata(h2copath+'APEX_H2CO_merge_high_smooth_noise.fits')
    f[0].data /= n
    f.writeto(h2copath+'APEX_H2CO_303_202_peaksn_smooth.fits',clobber=True)

def docleannhits():
    """ not really used now """
    f = fits.open(h2copath+'APEX_H2CO_merge_high_nhits.fits')
    nh = f[0].data
    nhm = scipy.ndimage.median_filter(nh, 5)
    f[0].data = nhm

def ph2cogrid(ntemp=50, trange=[10,200], abundances=(10**-8.5,10**-9),
              nh2=3e22, densityrange=[4,7]):
    import pyradex

    temperatures=np.linspace(trange[0],trange[1],ntemp)

    # initial density; will be modified later
    density = 1e4

    deltav = 5.0 # km/s

    R = pyradex.Radex(species='ph2co-h2',
                      abundance=abundances[0],
                      collider_densities={'H2':density},
                      deltav=deltav,
                      column=None,
                      temperature=temperatures[0],
                      h2column=nh2)

    Xarr = {}
    for abundance in abundances:
        Xarr[abundance] = {}
        for nh2 in (nh2,):

            densities = [10**x for x in xrange(*densityrange)]
            ratio1 = {d:[] for d in densities}
            ratio2 = {d:[] for d in densities}
            f1 = {d:[] for d in densities}
            f2 = {d:[] for d in densities}
            f3 = {d:[] for d in densities}

            for density in densities:
                R.density = {'H2': density}
                for temperature in temperatures:
                    R.temperature = temperature
                    print R.run_radex(),

                    F1 = R.T_B[2]  # 218.222192 3_0_3
                    F2 = R.T_B[12] # 218.760066 3_2_1
                    F3 = R.T_B[9]  # 218.475632 3_2_2

                    ratio1[density].append(F2/F1)
                    ratio2[density].append(F3/F1)
                    f3[density].append(F3)
                    f2[density].append(F2)
                    f1[density].append(F1)
                print

            f1 = {d:np.array([x.value for x in f1[d]]) for d in densities}
            f2 = {d:np.array([x.value for x in f2[d]]) for d in densities}
            f3 = {d:np.array([x.value for x in f3[d]]) for d in densities}
            ratio1 = {d:np.array(ratio1[d]) for d in densities}
            ratio2 = {d:np.array(ratio2[d]) for d in densities}

            Xarr[abundance][nh2] = {'flux1':f1,
                                    'flux2':f2,
                                    'flux3':f3,
                                    'ratio1':ratio1,
                                    'ratio2':ratio2}

    return Xarr

class TemperatureMapper(object):
    """
    For lazier evaluation of temperature mapping function
    """
    def __init__(self, trange=[10,300], ntemp=100, **kwargs):
        self.trange = trange
        self.ntemp = ntemp
        self.kwargs = kwargs

    def init(self):
        self.Xarr = ph2cogrid(trange=self.trange, ntemp=self.ntemp,
                              abundances=(10**-8.5,), nh2=3e22, **self.kwargs)
        self.temperatures = np.linspace(self.trange[0], self.trange[1], self.ntemp)


    def get_mapper(self, lineid, tmin=np.nan, tmax=np.nan):
        if not hasattr(self,'temperatures'):
            self.init()

        rationame = {'321220': 'ratio1',
                     '322221': 'ratio2'}[lineid]

        # ugly hack because ph2co is indexed with floats
        # Use FIXED abundance, FIXED column, FIXED density
        ratios = self.Xarr[self.Xarr.keys()[0]][3e22][rationame][1e4]

        def ratio_to_tem(r):
            inds = np.argsort(ratios)
            return np.interp(r, ratios[inds], self.temperatures[inds], tmin,
                             tmax)

        return ratio_to_tem

    def __call__(self, x, lineid='321220', **kwargs):
        return self.get_mapper(lineid, **kwargs)(x)

if 'tm' not in locals():
    tm = TemperatureMapper(trange=[10,300],ntemp=100)

def do_temperature(ratio=True):
    temperaturemap(tm, ratio=ratio)

def temperaturemap(ratio_to_tem, path=h2copath, ratio=True):
    if ratio:
        doratio()

    import scipy.stats

    for suf_ in ('{0}','{0}_integ','_cube{0}'):
        for smooth in ('_bl','_smooth_bl','_vsmooth_bl'):

            suf = suf_.format(smooth)

            for highline in ('321220','322221'):
                pfx = '{2}/H2CO_{0}_to_303202{1}'.format(highline,suf,path)
                if os.path.exists(pfx+'.fits'):
                    log.info("Temperature mapping {0}".format(pfx))
                else:
                    log.info("Skipping {0}".format(pfx))
                    continue
                rmap = fits.getdata(pfx+'.fits')
                tmap = ratio_to_tem(rmap, highline)#, tmin=10, tmax=300)
                rf = fits.open(pfx+'.fits')
                rf[0].header['BUNIT'] = 'K'
                rf[0].header['BTYPE'] = 'TKIN'
                del rf[0].header['LONPOLE']
                del rf[0].header['LATPOLE']
                rf[0].data = tmap
                rf.writeto(pfx+'_temperature.fits',clobber=True)
                #mask = fits.getdata('APEX_H2CO_303_202_smooth_mask_integ.fits') > 0.018
                if 'cube' not in suf:
                    mask = fits.getdata(path+'APEX_H2CO_322_221{0}_mask_integ.fits'.format(smooth)) > 0.0025
                else:
                    mask = fits.getdata(path+'APEX_H2CO_303_202{0}_mask.fits'.format(smooth.rstrip('_bl'))).astype('bool')
                tmap[True-mask] = np.nan
                rf[0].data = tmap
                rf.writeto(pfx+'_temperature_masked.fits',clobber=True)

                if 'cube' in suf:
                    rf[0].data = np.nanmax(tmap, axis=0)
                    rf.writeto(pfx+'_peaktemperature.fits', clobber=True)
                    rf[0].data = scipy.stats.nanmedian(tmap, axis=0)
                    rf.writeto(pfx+'_midtemperature.fits', clobber=True)

                    wt = fits.getdata(path+'APEX_H2CO_303_202{0}.fits'.format(smooth)).astype('bool')
                    wt[(True-mask) | (tmap==0) | (True-np.isfinite(tmap))] = 0
                    rf[0].data = np.nansum(tmap*wt, axis=0)/np.nansum(wt, axis=0)
                    rf.writeto(pfx+'_wtdmeantemperature.fits', clobber=True)

def mask_out_ch3oh(smooth='_smooth', dpath=mergepath):
    nu_ch3oh = all_lines['CH3OH_422_312']
    nu_h2co = all_lines['H2CO_322_221']
    v_ch3oh = ((nu_ch3oh - nu_h2co)/nu_h2co * constants.c).to(u.km/u.s).value

    hdu = fits.open(dpath+'APEX_H2CO_322_221{0}.fits'.format(smooth))[0]
    dv = hdu.header['CDELT3']
    shift = v_ch3oh / dv
    log.info("CH3OH Masking: dv: {0} shift: {1} ".format(dv,shift))

    mask = fits.getdata(dpath+'APEX_H2CO_303_202{0}_mask.fits'.format(smooth)).astype('bool')
    log.info("Mask shape: {0}".format(mask.shape))
    newmask = mask*False
    log.info("NewMask shape: {0}".format(newmask.shape))
    newmask[np.abs(shift):,:,:] = mask[:-np.abs(shift),:,:]
    log.info("NewMask number of masked pixels: {0}".format(newmask.sum()))
    hdu.data[newmask] = np.nan
    hdu.writeto(dpath+'APEX_H2CO_322_221{0}_CH3OHchomped.fits'.format(smooth), clobber=True)

    hdu.data[True-mask] = np.nan
    hdu.writeto(dpath+'APEX_H2CO_322_221{0}_CH3OHchomped_masked.fits'.format(smooth), clobber=True)

    integrate_mask(dpath+'APEX_H2CO_322_221{0}_CH3OHchomped'.format(smooth),
                   mask=dpath+'APEX_H2CO_303_202{0}_mask.fits'.format(smooth))

def do_mask_ch3oh(dpath=mergepath):
    mask_out_ch3oh('', dpath=dpath)
    # spatial smoothing = 2pix
    mask_out_ch3oh('_smooth', dpath=dpath)
    # spatial smoothing = 4pix
    mask_out_ch3oh('_vsmooth', dpath=dpath)

    mask_out_ch3oh('_bl', dpath=dpath)
    # spatial smoothing = 2pix
    mask_out_ch3oh('_smooth_bl', dpath=dpath)
    # spatial smoothing = 4pix
    mask_out_ch3oh('_vsmooth_bl', dpath=dpath)

def do_2014(datasets=datasets_2014, scanblsub=False):
    #datasets = ['E-093.C-0144A.2014APR02/E-093.C-0144A-2014-2014-04-01',
    #            'E-093.C-0144A.2014APR03/E-093.C-0144A-2014-2014-04-02']
    #build_cube_2014('MAP_001', datasets=datasets, scanblsub=True, lowhigh='low')
    #build_cube_2014('MAP_001', datasets=datasets, scanblsub=True, lowhigh='high')
    #build_cube_2014('MAP_001', datasets=datasets, scanblsub=False, lowhigh='high_nosub')

    for dataset in datasets:
        for source in datasets[dataset]:
            build_cube_2014(source, datasets=[dataset], scanblsub=scanblsub,
                            outpath=mergepath,
                            datapath=april2014path,
                            lowhigh='low')
            build_cube_2014(source, datasets=[dataset], scanblsub=scanblsub,
                            outpath=mergepath,
                            datapath=april2014path,
                            lowhigh='high')


def do_2014_merge(datasets=datasets_2014):
    log.info("Starting merge")
    for lowhigh in ('high','low',):
        mergefile = 'APEX_H2CO_2014_merge_{0}'.format(lowhigh)
        log.info("Making blanks")
        lowest_freq = 218.4e9 if lowhigh=='high' else 216.9e9
        make_blanks_merge(os.path.join(mergepath,mergefile), lowhigh=lowhigh,
                          lowest_freq=lowest_freq, width=2.5*u.GHz)
        mapnames = ['MAP_{0:03d}'.format(ii) for ii in range(1,130)]
        log.info("Building cubes: "+str(mapnames)+" "+lowhigh)
        build_cube_2014(mapnames,
                        mergefile=mergefile,
                        outpath=mergepath,
                        datapath=april2014path,
                        lowhigh=lowhigh,
                        datasets=datasets)

def get_info_2014(datapath='/Users/adam/work/h2co/apex/april2014/',
                  datasets=datasets_2014):
    info = {}
    for dataset in datasets:
        spectra,headers,indices = load_apex_cube(apex_filename=os.path.join(datapath,dataset)+".apex",
                                                 skip_data=True)
        info[dataset] = set([h['OBJECT'] for h in headers])
        log.info("{0}:{1}".format(dataset, str(info[dataset])))

    return info


def identify_scans_fromcoords(gal):
    # identify where the *derivative* changes signs
    # each np.diff shifts 1 to the left
    # 2 np.diffs -> +2 to index
    scans = 2+np.where(np.diff(np.sign(np.diff(gal.l.wrap_at(180*u.deg)))))[0]
    return scans

def subtract_scan_linear_fit(data, scans, mask_pixels=None,
                             verbose=False, smoothing_width=10,
                             automask=False, smooth_all=False,
                             smoothing_kernel_size_scale=40,
                             nsigma_ignore=1, return_mask=False):
    """
    Use linear algebra to fit a time-baseline to each scan to remove spectral
    baseline drifts.

    WARNING: This may remove map-spanning signals!!  That can be BAD for 13CO!

    Source:
    http://stackoverflow.com/questions/20343500/efficient-1d-linear-regression-for-each-element-of-3d-numpy-array
    (includes a solution for masked arrays: this will be EXTREMELY useful!)

    Parameters
    ----------
    data : np.ndarray
        2D data, with time along axis 0 and frequency along axis 1
    scans : np.ndarray
        The endpoints of the scans.  Should not include 0 or naxis
    divscale : bool
        DISABLED: this is crazy
        If True, will use only the slope and will divide out the normalized
        slope rather than subtracting
    mask_pixels : None or np.ndarray
        A mask array to select pixels to interpolate the fits across in
        the *Frequency* axis
    automask : bool
        Mask any scans with a mean > the overall mean + 1 stddev.  The data are
        slightly smoothed first if automask > 1.
    verbose : bool
        Print out simple stats about the fits
    smoothing_kernel_size_scale : int
        The size multiplier of the smoothing kernel used for interpolation in
        the frequency domain; smoothing_kernel_size_scale * smoothing_width
        defines the number of pixels to use when interpolating
    nsigma_ignore : float
        Fit masking control parameter.  Pixels with values greater than the
        mean noise + nsigma_ignore * std(mean_spectrum) will be ignored for
        fitting then interpolated back over later
    return_mask : bool
        Return an array of the mask used for each scan
    """
    
    #dmeans = data[:,percentile*data.shape[1]:(1-percentile)*data.shape[1]].mean(axis=1)

    dsub = data*0

    timeaxis = 0
    freqaxis = 1

    # Kernel must be ODD
    kernel_size = smoothing_kernel_size_scale * smoothing_width
    if kernel_size % 2 == 0:
        kernel_size += 1

    masklist = []

    for ii,jj in zip([0]+scans.tolist(),
                     scans.tolist()+[data.shape[timeaxis]]):
        x = np.arange(jj-ii)

        if automask:
            mean_spectrum = data[ii:jj,:].mean(axis=timeaxis)
            if automask > 1:
                mean_spectrum = convolve(mean_spectrum,
                                         Gaussian1DKernel(stddev=automask))
            mask_pixels = (mean_spectrum < (mean_spectrum.mean() +
                                            nsigma_ignore*mean_spectrum.std()))
            if verbose:
                nflag = (~mask_pixels).sum()
                log.info(("Masked {0} pixels for scanblsub fitting"
                          " in scan {1}-{2} "
                          "({3}%)").format(nflag, ii, jj,
                                           nflag/float(mask_pixels.size),)
                          )

        if mask_pixels is None:
            y = data[ii:jj,:]
        else:
            # mask_pixels is an include mask
            inds = np.arange(data.shape[freqaxis])[mask_pixels]
            y = data[ii:jj,mask_pixels]
            if return_mask and automask > 0:
                masklist.append(mask_pixels)

        # X is a vector of the X-values and a constant (1)
        # Becomes set of equations y = m x + b  ||  y = X mb
        X = np.c_[x,np.ones(jj-ii)]
        mb = np.linalg.lstsq(X,y)[0]
        
        if mask_pixels is not None:
            # Mask out the bad values, interpolate using a wide gaussian that
            # ignores nans
            m = np.zeros(data.shape[freqaxis]) + np.nan
            m[inds] = mb[0,:]
            m = convolve(m, Gaussian1DKernel(stddev=smoothing_width,
                                             x_size=kernel_size))

            b = np.zeros(data.shape[freqaxis]) + np.nan
            b[inds] = mb[1,:]
            b = convolve(b, Gaussian1DKernel(stddev=smoothing_width,
                                             x_size=kernel_size))

            # restore initial sampling unless we want smooth
            if not smooth_all:
                m[inds] = mb[0,:]
                b[inds] = mb[1,:]

            mb = np.array([m,b])

        dsub[ii:jj,:] = data[ii:jj,:] - np.inner(X,mb.T)

    log.info("Fit {0} scans with mean slopes {1} and offset {2}".format(len(scans)+1,
                                                                        mb.mean(axis=1)[0],
                                                                        mb.mean(axis=1)[1]))
    if np.any(np.isnan(dsub)):
        warnings.warn("There were NaNs left over from time-baseline subtraction.")
        dsub[np.isnan(dsub)] = 0

    if return_mask:
        return dsub, np.array(masklist)

    return dsub

def efuncs(arr, neig=None, return_others=False, huge_limit=500):
    """
    Determine eigenfunctions of an array for use with
    PCA cleaning

    Parameters
    ----------
    arr : `numpy.ndarray`
        The array (2D)
    neig : None or int
        The number of eigenvalues to compute.  Smaller = faster!
        None = All!
    huge_limit : int
        The limit above which an error will be raised (for large arrays, this
        can take *forever*)
    return_others : bool
        Return the evals, evects, and covmat or just the efuncs?

    Returns
    -------
    efuncarr : np.ndarray
        The eigenfunctions

    Optional Returns
    ----------------
    covmat : np.ndarray
        Symmetric covariance matrix
    evals : np.ndarray
        1D array of eigenvalues
    evects : np.ndarray
        Eigenvectors
    """
    if hasattr(arr,'filled'):
        arr = arr.filled(0)
    if arr.shape[1] > huge_limit and not neig:
        log.critical("Very large eigenvalue computation!"
                     " Danger stranger! Stranger danger!")
        import ipdb; ipdb.set_trace()
    covmat = np.dot(arr.T.conj(),arr)

    # assert covariance matrix is Hermitian
    # (symmetric under transpose + conjugation)
    if not (covmat.T.conj() == covmat).all():
        diff = (covmat.T.conj() - covmat)
        worst_diff_ind = np.argmax(np.abs(diff))
        worst_diff = diff.flat[worst_diff_ind]/covmat.flat[worst_diff_ind]
        log.warning("There are differences between the upper "
                    "and lower triangular components of the "
                    "covariance matrix; this is probably a "
                    "floating point error and should not be terrible."
                    "  The relative error is {wd}.".format(wd=worst_diff))
        if np.abs(worst_diff) > 1e-4:
            log.warning("Actually, that's a pretty large error.  "
                        "You may be in trouble.")
    
    # Changed from np.linalg.eig to scipy.linalg.eigh
    # and numpy.linalg.eigh, which both return values in
    # the opposite order from np.linalg.eig
    if neig:
        sz = covmat.shape[1]
        eva, eve = scipy.linalg.eigh(covmat,
                                     eigvals=(sz-neig,sz-1))
        # eigh returns values in opposit order from np.linalg.eig
        # we also want a fully populated matrix so the size stays
        # the same
        inds = np.argsort(eva)[::-1]

        evals = np.zeros(sz)
        evals[:neig] = eva[inds]
        evects = np.zeros([sz,sz])
        evects[:, :neig] = eve[:,inds]

    else:
        evals,evects = np.linalg.eigh(covmat)
        inds = np.argsort(evals)[::-1]
        evals = evals[inds]
        evects = evects[:,inds]

    efuncarr = np.dot(arr,evects)
    if return_others:
        return efuncarr,covmat,evals,evects
    else:
        return efuncarr

def PCA_clean(data,
              smoothing_scale=25., # should be ~200 for SEDIGISM
              timeaxis=0,
              freqaxis=1,
              ncomponents=3,
              diagplotfilename=None,
              scans=None,
              maxntimes=5000,
             ):
    """
    Remove N PCA components in the time direction

    TODO: speed up by downsampling in TIME as well; we don't expect large
    second-to-second variations REVISE: No, actually, there are sharp
    jumps in time.

    Maybe scan-by-scan pca is faster?

    Smoothing scale is ~200 in total, which means 25 for pre-downsampled
    CMZ data

    Parameters
    ----------
    data : `numpy.ndarray`
        2D data, with dimensions ``[times, frequencies]`` (or reversed if
        ``timeaxis`` and ``freqaxis`` are appropriately specified)
    smoothing_scale : float
        The scale over which frequencies should be smoothed prior to performing
        the PCA analysis.  This is the width of a gaussian.  The data will be
        downsampled by a factor (1/5)*smoothing_scale
    timeaxis : int
    freqaxis : int
        The axis #'s of the frequency and time data
    ncomponents : int
        The number of PCA components to remove.  3 is empirically decent, but
        it's very important to test this #
    diagplotfilename : None or str
        A filename to save a diagnostic plot in.  The plot shows the first
        ``ncomponents`` eigenfunctions.
    scans : list
        A list of scans.  If these are specified, the PCA analysis will be done
        on a scan-by-scan basis, in which the most-correlated N components will
        be identified in each scan.  This is not obviously the best thing to
        do, but it can be useful.
    maxntimes : int or None
        If specified, the timestream will be chunked out into sections with
        length < maxntimes before doing PCA computations.  In principle, this
        can be used to overcome memory limitations, but it should be used with
        caution as the locations of the splits are somewhat arbitrary and could
        result in different principle component selections if the data aren't
        well-behaved.
    """

    if freqaxis == 0 and timeaxis == 1:
        data = data.swapaxes(0,1)
    elif freqaxis != 1 or timeaxis != 0:
        raise ValueError("Invalid axis specification.")

    if np.any(np.isnan(data)):
        warnings.warn("There were NaNs in the PCA-target data")
        import ipdb; ipdb.set_trace()
        data = np.nan_to_num(data)

    if maxntimes and scans is None:
        ntimes = data.shape[0]
        if ntimes > maxntimes:
            nsplits = np.ceil(ntimes/float(maxntimes))
            length = ntimes/nsplits
            # Split with equal length, but leave out the starting point
            # and the end point since those are both added
            splits = np.linspace(0, ntimes, nsplits+1)[1:-1]
            scans = splits.astype('int')

    if scans is not None:
        all_data = data
        all_dsub = np.empty(data.shape)
        for start,end in zip([0]+scans.tolist(),
                             scans.tolist()+[data.shape[0]]):
            log.info("Computing PCA on an array with shape"
                     " {0}".format(data[start:end,:].shape))
            dsub,efuncarr = PCA_subtract(data[start:end,:],
                                         smoothing_scale=smoothing_scale,
                                         ncomponents=ncomponents)
            if start == 0:
                efuncs = efuncarr[:,:ncomponents]
            else:
                efuncs += efuncarr[:,:ncomponents]
            all_dsub[start:end,:] = dsub
        dsub = all_dsub
        efuncarr = efuncs / (len(scans)+1.) # Average removed efuncs
    else:
        log.info("Computing PCA on an array with shape"
                 " {0}".format(data.shape))
        dsub,efuncarr = PCA_subtract(data,
                                     smoothing_scale=smoothing_scale,
                                     ncomponents=ncomponents)


    if diagplotfilename is not None:
        fig = pl.figure(4)
        fig.clf()
        ax = fig.gca()
        for ii in range(ncomponents):
            ax.plot(efuncarr[:,ii], label=str(ii), linewidth=2, alpha=0.5)
        ax.legend(loc='best')
        checkdir_makedir(diagplotfilename)
        fig.savefig(diagplotfilename, bbox_inches='tight')

    if freqaxis == 0 and timeaxis == 1:
        dsub = dsub.swapaxes(0,1)

    return dsub.real

def PCA_subtract(data, smoothing_scale=None, ncomponents=3):
    """

    Parameters
    ----------
    data : `numpy.ndarray`
        2D data, with dimensions (times, frequencies)
    smoothing_scale : float
        The scale over which frequencies should be smoothed prior to performing
        the PCA analysis.  This is the width of a gaussian.  The data will be
        downsampled by a factor (1/5)*smoothing_scale

    Returns
    -------
    dsub : `numpy.ndarray`
        The data with ``ncomponents`` principle components removed
    efuncarr :
    """
    t0 = time.time()
    log.info("PCA will remove {0} components".format(ncomponents))
    if smoothing_scale:
        log.info(("PCA cleaning an image with size {0},"
                  " which will downsample to {1}").format(data.shape,
                                                          (data.shape[0],
                                                           data.shape[1]/(smoothing_scale/5))))

        sm_data = filters.gaussian_filter1d(data, smoothing_scale,
                                            axis=1, mode='mirror').real

        efuncarr,covmat,evals,evects = efuncs(sm_data[:,::smoothing_scale/5].T,
                                              neig=ncomponents,
                                              huge_limit=1000,
                                              return_others=True)
    else:
        log.info("PCA cleaning an image with size {0}".format(data.shape))
        
        efuncarr,covmat,evals,evects = efuncs(data.T,
                                              neig=ncomponents,
                                              huge_limit=1000,
                                              return_others=True)

    log.info("Completed PCA (eigenfunction/vector) computation"
             " in {0} seconds.".format(time.time()-t0))

    # Zero-out the components we want to keep
    # (technically no longer necessary: this should be a null operation)
    efuncarr[:,ncomponents:] = 0

    to_subtract = np.inner(efuncarr,evects).T

    if smoothing_scale:
        ifunc = interpolate.interp1d(np.arange(to_subtract.shape[1]),
                                     to_subtract,
                                     axis=1)
        to_subtract = ifunc(np.linspace(0, to_subtract.shape[1]-1, data.shape[1]))

    dsub = data - to_subtract

    return dsub, efuncarr

def extract_co_subcubes(mergepath=april2014path):
    extract_subcube(os.path.join(mergepath,'APEX_H2CO_2014_merge_high.fits'),
                    os.path.join(mergepath,'APEX_13CO_2014_merge.fits'),
                    linefreq=220.39868*u.GHz, naxis3=500, vmin=-225*u.km/u.s,
                    vmax=275*u.km/u.s)
    extract_subcube(os.path.join(mergepath,'APEX_H2CO_2014_merge_high.fits'),
                    os.path.join(mergepath,'APEX_C18O_2014_merge.fits'),
                    linefreq=219.56036*u.GHz, naxis3=500, vmin=-225*u.km/u.s,
                    vmax=275*u.km/u.s)
