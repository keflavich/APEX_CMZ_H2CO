import pyspeckit
from pyspeckit.spectrum.readers import read_class
from pyspeckit import cubes
import numpy as np
from astropy import wcs
from astropy import coordinates
from astropy import units as u
from astropy import constants
from astropy.utils.console import ProgressBar
from astropy.convolution import convolve, Gaussian1DKernel
from sdpy import makecube
from astropy.io import fits
from FITS_tools import cube_regrid
from FITS_tools.load_header import get_cd
from astropy.wcs import WCS
import FITS_tools
import scipy.ndimage
import time
import mpl_plot_templates
import pylab as pl
import os
from astropy import log
import glob
from scipy.ndimage import filters
from scipy import signal,interpolate
import image_tools

datasets_2014 = {'E-093.C-0144A.2014JUN01/E-093.C-0144A-2014-2014-05-31': ('MAP_007',),
                 'E-093.C-0144A.2014MAY30/E-093.C-0144A-2014-2014-05-29': ('MAP_002','MAP_003','MAP_004'),
                 'E-093.C-0144A.2014MAY31/E-093.C-0144A-2014-2014-05-30': ('MAP_005','MAP_006'),
                 'E-093.C-0144A.2014APR02/E-093.C-0144A-2014-2014-04-01': ('MAP_001',),
                 'E-093.C-0144A.2014APR03/E-093.C-0144A-2014-2014-04-02': ('MAP_001',),
                 'E-093.C-0144A.2014JUN02/E-093.C-0144A-2014-2014-06-01': ('MAP_009','MAP_010','MAP_008',),
                 'E-093.C-0144A.2014JUN03/E-093.C-0144A-2014-2014-06-02': ('MAP_011','MAP_012','MAP_013', 'MAP_018', 'MAP_019'),
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
h2copath = '/Users/adam/work/h2co/apex/h2co_cubes/'
mergepath = '/Users/adam/work/h2co/apex/merged_datasets/'
aorawpath = '/Users/adam/work/h2co/apex/2010_reduced/2010_raw/'
aopath = '/Users/adam/work/h2co/apex/2010_reduced/'
diagplotdir = '/Users/adam/work/h2co/apex/diagnostic_plots/'

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

    spectra,headers,indices,data,hdrs,gal = load_2013_dataset_for_debugging(skip_data=False, lowhigh='high')
    make_blanks_freq(gal, hdrs[0], test, clobber=True)

    dmeansub,gal,hdrs = process_data(data, gal, hdrs, dataset=test,
                                     subspectralmeans=True, scanblsub=False)
    add_apex_data(dmeansub, hdrs, gal, test, retfreq=True, varweight=True,)
    make_blanks_freq(gal, hdrs[0], test+"_blsub", clobber=True)
    dscube = cube_regrid.downsample_cube(fits.open(test+".fits")[0], factor=4)
    dscube.writeto(test+"_ds.fits",clobber=True)

    dspecsub,gal,hdrs = process_data(data, gal, hdrs, dataset=test+"_blsub",
                                     subspectralmeans=True, scanblsub=True)
    add_apex_data(dspecsub, hdrs, gal, test+"_blsub", retfreq=True, varweight=True,)
    dscube = cube_regrid.downsample_cube(fits.open(test+"_blsub.fits")[0], factor=4)
    dscube.writeto(test+"_blsub_ds.fits",clobber=True)

    freq = hdr_to_freq(hdrs[0])
    mask = make_line_mask(freq)

    return spectra,headers,indices,data,hdrs,gal,dspecsub,dmeansub,freq,mask

def load_2013_dataset_for_debugging(lowhigh='low', downsample_factor=8,
                                    dataset='M-091.F-0019-2013-2013-06-11',
                                    datapath=june2013datapath,
                                    xscan=37986,
                                    skip_data=True):
    """
    Example:

    spectra,headers,indices, data,hdrs,gal = load_2013_dataset_for_debugging(skip_data=False)
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
    xtel = 'AP-H201-X202' if lowhigh=='low' else 'AP-H201-X201'

    apex_filename=datapath+dataset+".apex"

    spectra,headers,indices = load_apex_cube(apex_filename,
                                             skip_data=skip_data,
                                             downsample_factor=downsample_factor)
    data, hdrs, gal = select_apex_data(spectra, headers, indices,
                                       sourcename='SGRA',
                                       shapeselect=4096,
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

def process_data(data, gal, hdrs, dataset, scanblsub=True,
                 subspectralmeans=True, verbose=False, noisefactor=1.5,
                 linemask=False, automask=2,
                 zero_edge_pixels=0,
                 pca_clean=False,
                 pcakwargs={},
                 **kwargs):

    timeaxis = 0
    freqaxis = 1

    log.info("Processing {0}".format(dataset))

    if zero_edge_pixels:
        # Force the Nth first/last frequency pixels to zero
        data[:,:zero_edge_pixels] = 0
        data[:,-zero_edge_pixels:] = 0

    if subspectralmeans:
        data -= data.mean(axis=freqaxis)[:,None]

    obsids = np.array([h['XSCAN'] for h in hdrs])

    if scanblsub:

        data_diagplot(data, dataset+"_presub",
                      **kwargs)
        for ii,xscan in enumerate(np.unique(obsids)):
            match = obsids == xscan
            # maybe mask=mask_pix.max(axis=timeaxis), ?
            #mask=mask_pix[ii], 
            data_diagplot(data[match], dataset+"_presub_obs%i" % xscan,
                          **kwargs)

        scans = identify_scans_fromcoords(gal)
        freq = hdr_to_freq(hdrs[0])
        if linemask:
            mask = make_line_mask(freq)
        else:
            mask = None
        dsub,mask_pix = subtract_scan_linear_fit(data, scans, mask_pixels=mask,
                                                 verbose=verbose,
                                                 automask=automask,
                                                 smooth_all=True,
                                                 return_mask=True)
        mask = mask_pix.max(axis=timeaxis).astype('bool')
    else:
        freq = None
        mask = None
        dsub = data

    if pca_clean:
        dsub = PCA_clean(dsub, **pcakwargs)

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
             **kwargs)

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

    diagplot(dsub, tsys, noise, dataset, freq=freq, mask=mask, **kwargs)
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

    makecube.add_data_to_cube(cubefilename+".fits", data=data, flatheader='header.txt',
                              cubeheader='cubeheader.txt', linefreq=218.22219, allow_smooth=True,
                              nhits=nhits,
                              data_iterator=data_iterator,
                              coord_iterator=coord_iterator,
                              velo_iterator=velo_iterator, debug=debug,
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

    makecube.generate_header(np.mean(lrange), np.mean(brange), naxis1=naxis1,
                             naxis2=naxis2, naxis3=4096, coordsys='galactic',
                             ctype3='VRAD',
                             bmaj=bmaj.to(u.deg).value,
                             bmin=bmaj.to(u.deg).value,
                             pixsize=pixsize.to(u.arcsec).value,
                             cunit3='km/s',
                             output_flatheader='header.txt',
                             output_cubeheader='cubeheader.txt',
                             cd3=header['VRES'],
                             crval3=-1*header['VRES']*header['RCHAN'],
                             crpix3=1,
                             clobber=True, bunit="K",
                             restfreq=restfreq.to(u.Hz).value, radio=True)

    makecube.make_blank_images(cubefilename, clobber=clobber)

def make_blanks_freq(gal, header, cubefilename, clobber=True, pixsize=7.2*u.arcsec):
    """ complete freq covg """

    lrange = gal.l.wrap_at(180*u.deg).deg.min()+15/3600.,gal.l.wrap_at(180*u.deg).deg.max()+15/3600.
    brange = gal.b.deg.min()+15/3600.,gal.b.deg.max()+15/3600.
    print "Map extent: %0.2f < l < %0.2f,  %0.2f < b < %0.2f" % (lrange[0],
                                                                 lrange[1],
                                                                 brange[0],
                                                                 brange[1])

    naxis1 = int((lrange[1]-lrange[0])/(pixsize.to(u.deg).value)+10)
    naxis2 = int((brange[1]-brange[0])/(pixsize.to(u.deg).value)+10)
    restfreq = (header['RESTF']*u.MHz)
    bmaj = (1.22*10*u.m / restfreq.to(u.m,u.spectral()))*u.radian
    rchan = header['RCHAN']

    #scalefactor = 1./downsample_factor
    #crpix3 = (rchan-1)*scalefactor+0.5+scalefactor/2.

    makecube.generate_header(np.mean(lrange), np.mean(brange), naxis1=naxis1,
                             naxis2=naxis2, naxis3=header['NCHAN'],
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
                             restfreq=restfreq.to(u.Hz).value, radio=True)

    makecube.make_blank_images(cubefilename, clobber=clobber)


def make_blanks_merge(cubefilename, lowhigh='low', clobber=True,
                      width=1.0*u.GHz, lowest_freq=None, pixsize =
                      7.2*u.arcsec):
    # total size is 2.3 x 0.4 degrees
    # 1150x
    # center is 0.55 -0.075
    naxis1 = 1150
    naxis2 = 200
    restfreq = 218222.192*u.MHz
    bmaj = (1.22*10*u.m / restfreq.to(u.m,u.spectral()))**-1*u.radian
    cd3 = ((1*u.km/u.s)/constants.c * 218.2*u.GHz).to(u.Hz).value
    naxis3 = int(np.ceil(((width / (218.2*u.GHz) * constants.c) / (u.km/u.s)).decompose().value))
    if lowest_freq is None:
        lowest_freq = 216.8e9 if lowhigh=='low' else 218e9

    makecube.generate_header(0.55, -0.075, naxis1=naxis1,
                             naxis2=naxis2, naxis3=naxis3, coordsys='galactic',
                             bmaj=bmaj.to(u.deg).value,
                             bmin=bmaj.to(u.deg).value,
                             pixsize=pixsize.to(u.arcsec).value,
                             cunit3='Hz',
                             ctype3='FREQ',
                             output_flatheader='header.txt',
                             output_cubeheader='cubeheader.txt',
                             cd3=cd3,
                             crval3=lowest_freq,
                             crpix3=1,
                             clobber=True, bunit="K",
                             restfreq=restfreq.to(u.Hz).value, radio=True)

    makecube.make_blank_images(cubefilename, clobber=clobber)

def data_diagplot(data, dataset, ext='png', newfig=False,
                  max_size=1024, freq=None):
    log.info("Doing diagnostics in "+dataset)
    if newfig:
        pl.figure()
    else:
        pl.figure(1)
        pl.clf()
    if (np.isnan(data)).all():
        print "ALL data is NaN in ", dataset
        import ipdb; ipdb.set_trace()

    if np.any([d > max_size for d in data.shape]):
        # downsample to *not less than* max_size
        factors = [int(np.floor(d / max_size)) for d in data.shape]
        data = image_tools.downsample(data, min(factors))

    axis = mpl_plot_templates.imdiagnostics(data)

    if freq is not None:
        axis.set_xticklabels(np.interp(axis.get_xticks(),
                                       np.arange(freq.size),
                                       freq))
        axis.set_xlabel("Frequency")
    else:
        axis.set_xlabel("Channel #")
    
    axis.set_ylabel("Integration #")

    try:
        pl.savefig(os.path.join(diagplotdir, dataset+"_diagnostics."+ext),bbox_inches='tight')
    except Exception as ex:
        log.info(ex)
        print ex
    return axis

def diagplot(data, tsys, noise, dataset, freq=None, mask=None, ext='png',
             newfig=False):
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
    pl.savefig(os.path.join(diagplotdir, dataset+"_tsys."+ext),bbox_inches='tight')

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
    pl.savefig(os.path.join(diagplotdir, dataset+"_masked."+ext),bbox_inches='tight')

    data_diagplot(data, dataset, ext=ext, newfig=newfig, freq=freq)

def build_cube_generic(window, freq=True, mergefile=None, datapath='./',
                       outpath='./', datasets=[], scanblsub=True,
                       shapeselect=None,
                       sourcename=None,
                       tsysrange=[100,250],
                       excludefitrange=None,
                       downsample_factor=None,
                       pixsize=7.2*u.arcsec,
                       kernel_fwhm=10/3600.,
                       pca_clean=False,
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
                  datapath=aorawpath,
                  outpath=aopath,
                  datasets=['O-085.F-9311A-2010','E-085.B-0964A-2010'],
                  scanblsub=True,
                  verbose=False,
                  debug=False,
                  **kwargs):
    """
    TODO: comment!

    kwargs are passed to process_data
    """
    if window not in ('low','high'):
        raise ValueError()
    if mergefile:
        cubefilename=os.path.join(outpath,'APEX_H2CO_merge_%s' % window)
    elif freq:
        cubefilename=os.path.join(outpath,'APEX_H2CO_Ao_Freq_%s' % window)
    else:
        cubefilename=os.path.join(outpath,'APEX_H2CO_Ao_%s' % window)

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

    log.info("Data has been collected and flagged, now adding to cube.")

    for dataset in all_data:
        data = all_data[dataset]
        hdrs = all_hdrs[dataset]
        gal  = all_gal[dataset]

        data, gal, hdrs = process_data(data, gal, hdrs, dataset+"_"+xtel,
                                       scanblsub=scanblsub, verbose=verbose,
                                       **kwargs)

        add_apex_data(data, hdrs, gal, cubefilename,
                      excludefitrange=excludefitrange,
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
                    datasets=['M-091.F-0019-2013-2013-06-08',
                              'M-091.F-0019-2013-2013-06-11',
                              'M-091.F-0019-2013-2013-06-12',
                              'M-091.F-0019-2013-2013-06-13'],
                    scanblsub=True,
                    verbose=True, **kwargs):
    if mergefile:
        cubefilename=os.path.join(outpath,mergefile)
    else:
        cubefilename=os.path.join(outpath,
                                  'APEX_H2CO_2013_%s' % lowhigh)


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

    # need two loops to avoid loading too much stuff into memory
    for dataset in datasets:

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
                                       **kwargs)

        add_apex_data(data, hdrs, gal, cubefilename, retfreq=True,
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

def make_high_mergecube(datasets_2014=datasets_2014):
    mergefile2 = 'APEX_H2CO_merge_high'
    make_blanks_merge(os.path.join(mergepath,mergefile2), lowhigh='high')
    build_cube_ao(window='high', mergefile=True, freq=True, outpath=mergepath)
    build_cube_2013(mergefile=mergefile2,
                    outpath=mergepath,
                    lowhigh='high',
                    scanblsub=True)

    log.info("Starting merge")
    for lowhigh in ('low','high',):
        mapnames = ['MAP_{0:03d}'.format(ii) for ii in range(1,130)]
        log.info("Building cubes: "+str(mapnames)+" "+lowhigh)
        build_cube_2014(mapnames,
                        mergefile=mergefile2,
                        outpath=mergepath,
                        lowhigh=lowhigh,
                        datasets=datasets_2014)


def make_low_mergecube(datasets_2014=datasets_2014):
    mergefile1 = 'APEX_H2CO_merge_low'
    make_blanks_merge(os.path.join(mergepath,mergefile1), lowhigh='low')
    build_cube_ao(window='low', mergefile=True, freq=True, outpath=mergepath)
    build_cube_2013(mergefile=mergefile1,
                    outpath=mergepath,
                    lowhigh='low',
                    scanblsub=True)

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
    if isinstance(mask,str):
        mask = fits.getdata(mask)
    ffile = fits.open(prefix+'.fits')
    ffile[0].data *= mask

    integ1,hdr = cubes.integ(ffile, [0,ffile[0].shape[0]], average=np.nanmean)
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
        #mask._wcs = subcube1.wcs
        subcube = subcube1.with_mask(subcube1>0.05)#.with_mask(mask)
        if subcube.shape[0] == 1:
            # implies out of range
            continue
        mom0 = subcube.moment0()
        mom1 = subcube.moment1()
        mom2 = subcube.moment2()
        outfn = 'projections/'+filename.replace(".fits","_{line}_{mom}.fits")
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

def moment_mask_cube(prefix, noise=None, kernelsize=[2,2,2], grow=1,
                     sigmacut=3):
    ffile = fits.open(prefix+'.fits')
    cube = ffile[0].data
    if noise is None:
        noise = fits.getdata(prefix+'_noise.fits')

    t0 = time.time()
    print "Initiating cube smooth."
    smcube = cube_regrid.gsmooth_cube(cube, kernelsize, use_fft=False,
                                      kernelsize_mult=3)
    print "Completed cube smooth in %i seconds" % (time.time()-t0)
    mask = smcube > noise*sigmacut

    mask_grow = scipy.ndimage.morphology.binary_dilation(mask, iterations=1)

    ffile[0].data[True-mask_grow] = np.nan
    ffile[0].writeto(prefix+"_moment_masked.fits", clobber=True)

    ffile[0].data = mask_grow.astype('int')
    ffile[0].writeto(prefix+"_mask.fits", clobber=True)

def do_momcube_masking_hi(prefix=h2copath+'APEX_H2CO_303_202'):
    compute_noise_high(prefix)
    moment_mask_cube(prefix)
    integrate_slices_high(prefix+'_moment_masked')

def extract_subcube(cubefilename, outfilename, linefreq=218.22219*u.GHz,
                    debug=False, smooth=False, vsmooth=False, naxis3=400, crval3=50):
    t0 = time.time()
    log.info("Extracting subcube at {0} from {1} with smooth={2} and vsmooth={3}".format(linefreq, cubefilename, smooth, vsmooth))
    ffile = fits.open(cubefilename)
    # I don't know why this is necessary, but there were weird bugs showing up
    # if I did not add this (some parts of the cube defaulted to 3e-319)
    # Added June 16: float32->float64 because otherwise regridding can fail
    ffile[0].data = ffile[0].data.astype('float64')
    hdr = ffile[0].header

    cdk = 'CD3_3' if 'CD3_3' in hdr else 'CDELT3'

    if debug:
        xarr = (np.arange(hdr['NAXIS3'])+1-hdr['CRPIX3']) * hdr[cdk] + hdr['CRVAL3']
        log.debug("xarr min={0} max={1}".format(xarr.min(),xarr.max()))

    hdr['CTYPE3'] = 'VELO'
    cdfrq = hdr[cdk]
    crvfreq = hdr['CRVAL3']
    crpixf = hdr['CRPIX3']
    hdr[cdk] = -((cdfrq/crvfreq)*constants.c).to(u.km/u.s).value
    hdr['CRVAL3'] = 0.0
    hdr['CRPIX3'] = (linefreq.to(u.Hz).value-crvfreq)/cdfrq + crpixf
    hdr['CUNIT3'] = 'km/s'

    if debug:
        xarr = (np.arange(hdr['NAXIS3'])+1-hdr['CRPIX3']) * hdr[cdk] + hdr['CRVAL3']
        log.debug("xarr min={0} max={1}".format(xarr.min(),xarr.max()))

    for k in ['CTYPE3','CRVAL3',cdk,'CRPIX3','CUNIT3']:
        assert ffile[0].header[k] == hdr[k]

    outhdr = hdr.copy()
    outhdr[cdk] = 1.0
    outhdr['RESTFREQ'] = linefreq.to(u.Hz).value
    outhdr['RESTFRQ'] = linefreq.to(u.Hz).value
    outhdr['CRVAL3'] = crval3
    outhdr['CRPIX3'] = (naxis3-1)/2.
    outhdr['NAXIS3'] = naxis3

    if smooth:
        #cubesm = gsmooth_cube(ffile[0].data, [3,2,2], use_fft=True,
        #                      psf_pad=False, fft_pad=False)
        # smoothed with 2 pixels -> sigma=10", fwhm=23"
        # this is an "optimal smooth", boosting s/n and smoothing to 36"
        # resolution.
        kw = 2 if not vsmooth else 4
        cubesm = cube_regrid.spatial_smooth_cube(ffile[0].data, kw,
                                                 use_fft=False,
                                                 numcores=4)
        cubesm = cube_regrid.spectral_smooth_cube(cubesm, 3/2.35,
                                                  use_fft=False,
                                                  numcores=4)
        ffile[0].data = cubesm

        outhdr[cdk] = 3.0
        outhdr['CRVAL3'] = 50
        outhdr['CRPIX3'] = 70
        outhdr['NAXIS3'] = 140
    
    if debug:
        xarr = (np.arange(outhdr['NAXIS3'])+1-outhdr['CRPIX3']) * outhdr[cdk] + outhdr['CRVAL3']
        log.debug("xarr min={0} max={1}".format(xarr.min(),xarr.max()))
        xarr = (-xarr/3e5) * linefreq + linefreq
        log.debug("xarr min={0} max={1}".format(xarr.min(),xarr.max()))
        return hdr,outhdr

    newhdu = cube_regrid.regrid_cube_hdu(ffile[0], outhdr, order=1,
                                         prefilter=False)

    newhdu.writeto(outfilename, clobber=True)

    log.info("Completed cube extraction in {0} seconds.".format(time.time()-t0))
    
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



def make_line_mask(freqarr, lines=bright_lines):
    mask = np.ones(freqarr.size, dtype='bool')
    for ln,lf in lines.iteritems():
        bw = bandwidths[ln]
        wh = (lf*1e3-bw < freqarr) & (lf*1e3+bw > freqarr)
        mask[wh] = False
    return mask


def do_extract_subcubes():

    for line,freq in all_lines.iteritems():
        log.info("Extracting {0}".format(line))
        if freq < 218:
            cubefilename = os.path.join(mergepath,'APEX_H2CO_merge_low_sub.fits')
        else:
            cubefilename = os.path.join(mergepath,'APEX_H2CO_merge_high_sub.fits')

        header = fits.getheader(cubefilename)
        ww = wcs.WCS(header)
        wspec = ww.sub([wcs.WCSSUB_SPECTRAL])
        nax = header['NAXIS%i' % (ww.wcs.spec+1)]
        freqarr = wspec.wcs_pix2world(np.arange(nax),0)[0]
        if freq*1e9 > freqarr.min() and freq*1e9 < freqarr.max():
            extract_subcube(cubefilename, 'APEX_{0}_smooth.fits'.format(line),
                            linefreq=freq*u.GHz, smooth=True)
            extract_subcube(cubefilename, 'APEX_{0}_vsmooth.fits'.format(line),
                            linefreq=freq*u.GHz, smooth=True, vsmooth=True)
            extract_subcube(cubefilename, 'APEX_{0}.fits'.format(line),
                            linefreq=freq*u.GHz)
        else:
            log.info("Skipping line {0}".format(line))


def do_everything():
    make_high_mergecube()
    make_low_mergecube()
    os.chdir(mergepath)
    os.system('./APEX_H2CO_merge_high_starlink_custom.sh')
    os.chdir('../')
    do_extract_subcubes()
    compute_noise_high()
    compute_noise_high(mergepath+'APEX_H2CO_merge_high_smooth',[203,272])
    compute_noise_high(mergepath+'APEX_H2CO_merge_high_vsmoothds',[203,272])
    compute_noise_high(mergepath+'APEX_H2CO_303_202_vsmooth',[107,141])
    compute_noise_low()
    moment_mask_cube(mergepath+'APEX_H2CO_303_202',
                     noise=fits.getdata(mergepath+'APEX_H2CO_merge_high_sub_noise.fits'))
    moment_mask_cube(mergepath+'APEX_H2CO_303_202_smooth',
                     noise=fits.getdata(mergepath+'APEX_H2CO_merge_high_smooth_noise.fits'),
                     sigmacut=5)
    moment_mask_cube(mergepath+'APEX_H2CO_303_202_vsmooth',
                     noise=fits.getdata(mergepath+'APEX_H2CO_303_202_vsmooth_noise.fits'),
                     sigmacut=5)
    integrate_mask(mergepath+'APEX_H2CO_303_202_smooth',mask=mergepath+'APEX_H2CO_303_202_smooth_mask.fits')
    integrate_mask(mergepath+'APEX_H2CO_303_202_vsmooth',mask=mergepath+'APEX_H2CO_303_202_vsmooth_mask.fits')

    for fn in glob.glob(os.path.join(mergepath,'APEX_H2CO_30*fits')):
        try:
            os.symlink(fn,
                       os.path.join(h2copath,os.path.split(fn)[-1]))
        except OSError:
            log.debug("Skipped file {0} because it exists".format(fn))

    # Create H2CO masks?
    do_momcube_masking_hi()

    for line in all_lines:
        if os.path.exists(mergepath+'APEX_{0}.fits'.format(line)):
            integrate_mask(mergepath+'APEX_{0}'.format(line))
            integrate_mask(mergepath+'APEX_{0}_smooth'.format(line),
                           mask=mergepath+'APEX_H2CO_303_202_smooth_mask.fits')
            integrate_mask(mergepath+'APEX_{0}_vsmooth'.format(line),
                           mask=mergepath+'APEX_H2CO_303_202_vsmooth_mask.fits')
    for line in all_lines:
        if os.path.exists(mergepath+'APEX_{0}.fits'.format(line)):
            baseline_cube(mergepath+'APEX_{0}.fits'.format(line),
                          maskfn=mergepath+'APEX_H2CO_303_202_smooth_mask.fits')
            baseline_cube(mergepath+'APEX_{0}_smooth.fits'.format(line),
                          maskfn=mergepath+'APEX_H2CO_303_202_smooth_mask.fits')
            baseline_cube(mergepath+'APEX_{0}_vsmooth.fits'.format(line),
                          maskfn=mergepath+'APEX_H2CO_303_202_vsmooth_mask.fits')

    for fn in glob.glob(os.path.join(mergepath,'APEX_H2CO_32*fits')):
        try:
            os.symlink(fn,
                       os.path.join(h2copath,os.path.split(fn)[-1]))
        except OSError:
            log.debug("Skipped file {0} because it exists".format(fn))

    do_temperature()

def do_temperature():
    temperaturemap(tm)

def baseline_cube(cubefn, maskfn, order=5):
    from pyspeckit.cubes.cubes import baseline_cube
    f = fits.open(cubefn)
    cube = f[0].data
    mask = fits.getdata(maskfn).astype('bool')
    t0 = time.time()
    print "Baselining cube {0}...".format(cubefn),
    bc = baseline_cube(cube, order, cubemask=mask)
    print " done ({0} seconds)".format(time.time()-t0)
    f[0].data = bc
    f.writeto(cubefn.replace(".fits","_bl.fits"), clobber=True)


def do_everything_2013extrafreqs():
    build_cube_2013(lowhigh='low',
                    scanblsub=True)
    build_cube_2013(lowhigh='high',
                    scanblsub=True)
    #raise NotImplementedError
    #compute_noise_extras(lowhigh='low',pixrange=[0,4096])
    #compute_noise_extras(lowhigh='high',pixrange=[0,4096])


def doratio():
    """ I swapped top and bottom because that's what the models were set up for... """

    top = fits.getdata(h2copath+'APEX_H2CO_303_202_smooth_mask_integ.fits')
    bottom = fits.getdata(h2copath+'APEX_H2CO_322_221_smooth_CH3OHchomped_mask_integ.fits')
    
    ratio = bottom/top
    ratio[ratio<0.0] = np.nan
    ratio[ratio>0.64] = np.nan
    
    f = fits.open(h2copath+'APEX_H2CO_303_202_smooth_mask_integ.fits')
    
    f[0].data = ratio
    
    f.writeto(h2copath+'H2CO_322221_to_303202.fits',clobber=True)

    bottom = fits.getdata(h2copath+'APEX_H2CO_321_220_smooth_mask_integ.fits')
    ratio = bottom/top
    ratio[ratio<0.0] = np.nan
    ratio[ratio>0.64] = np.nan
    
    f[0].data = ratio
    
    f.writeto(h2copath+'H2CO_321220_to_303202.fits',clobber=True)

    ##### cube #####
    for smooth in ('smooth','vsmooth'):
        top = fits.getdata(h2copath+'APEX_H2CO_303_202_{0}.fits'.format(smooth))
        bottom = fits.getdata(h2copath+'APEX_H2CO_322_221_{0}_CH3OHchomped_masked.fits'.format(smooth))
        
        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>0.64] = np.nan
        
        f = fits.open(h2copath+'APEX_H2CO_303_202_{0}_mask.fits'.format(smooth))
        
        f[0].data = ratio
        
        f.writeto(h2copath+'H2CO_322221_to_303202_cube_{0}.fits'.format(smooth),clobber=True)

        bottom = fits.getdata(h2copath+'APEX_H2CO_321_220_{0}.fits'.format(smooth))
        
        ratio = bottom/top
        ratio[ratio<0.0] = np.nan
        ratio[ratio>0.64] = np.nan
        
        f = fits.open(h2copath+'APEX_H2CO_303_202_{0}.fits'.format(smooth))
        
        # mask out...
        f[0].data = ratio * f[0].data.astype('bool')
        
        f.writeto(h2copath+'H2CO_321220_to_303202_cube_{}.fits'.format(smooth),clobber=True)



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

def ph2cogrid(ntemp=50, trange=[10,200], abundance=10**-8.5,
              nh2=3e22):
    import pyradex

    temperatures=np.linspace(trange[0],trange[1],ntemp)

    # initial density; will be modified later
    density = 1e4

    deltav = 5.0 # km/s

    R = pyradex.Radex(species='ph2co-h2',
                      abundance=abundance,
                      collider_densities={'H2':density},
                      deltav=deltav,
                      column=None,
                      temperature=temperatures[0],
                      h2column=nh2)

    Xarr = {}
    for abundance in (abundance,):
        Xarr[abundance] = {}
        for nh2 in (nh2,):

            densities = [10**x for x in xrange(4,7)]
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

            Xarr[abundance][nh2] = (f1,f2,f3,ratio1,ratio2)

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
        self.Xarr = ph2cogrid(trange=self.trange, ntemp=self.ntemp, abundance=10**-8.5,
                              nh2=3e22, **self.kwargs)
        self.temperatures = np.linspace(self.trange[0], self.trange[1], self.ntemp)


    def get_mapper(self, tmin=np.nan, tmax=np.nan):
        if not hasattr(self,'temperature'):
            self.init()

        # ugly hack because ph2co is indexed with floats
        ratios = self.Xarr[self.Xarr.keys()[0]][3e22][3][1e4]

        def ratio_to_tem(r):
            inds = np.argsort(ratios)
            return np.interp(r, ratios[inds], self.temperatures[inds], tmin,
                             tmax)

        return ratio_to_tem

    def __call__(self,x, **kwargs):
        return self.get_mapper(**kwargs)(x)

if 'tm' not in locals():
    tm = TemperatureMapper(trange=[10,300],ntemp=100)

def temperaturemap(ratio_to_tem, path=h2copath):
    doratio()

    import scipy.stats

    for smooth in ('smooth','vsmooth'):
        for suf in ('','_cube_{0}'):
            if suf:
                suf = suf.format(smooth)
            for highline in ('321220','322221'):
                pfx = '{2}/H2CO_{0}_to_303202{1}'.format(highline,suf,path)
                rmap = fits.getdata(pfx+'.fits')
                tmap = ratio_to_tem(rmap)#, tmin=10, tmax=300)
                rf = fits.open(pfx+'.fits')
                rf[0].header['BUNIT'] = 'K'
                rf[0].header['BTYPE'] = 'TKIN'
                rf[0].data = tmap
                rf.writeto(pfx+'_temperature.fits',clobber=True)
                #mask = fits.getdata('APEX_H2CO_303_202_smooth_mask_integ.fits') > 0.018
                if not suf:
                    mask = fits.getdata('APEX_H2CO_322_221_{0}_mask_integ.fits'.format(smooth)) > 0.0025
                else:
                    mask = fits.getdata('APEX_H2CO_303_202_{0}_mask.fits'.format(smooth)).astype('bool')
                tmap[True-mask] = np.nan
                rf[0].data = tmap
                rf.writeto(pfx+'_temperature_masked.fits',clobber=True)

                if 'cube' in suf:
                    rf[0].data = np.nanmax(tmap, axis=0)
                    rf.writeto(pfx+'_peaktemperature.fits', clobber=True)
                    rf[0].data = scipy.stats.nanmedian(tmap, axis=0)
                    rf.writeto(pfx+'_midtemperature.fits', clobber=True)

                    wt = fits.getdata('APEX_H2CO_303_202_smooth.fits').astype('bool')
                    wt[(True-mask) | (tmap==0) | (True-np.isfinite(tmap))] = 0
                    rf[0].data = np.nansum(tmap*wt, axis=0)/np.nansum(wt, axis=0)
                    rf.writeto(pfx+'_wtdmeantemperature.fits', clobber=True)

def mask_out_ch3oh(smooth='smooth'):
    nu_ch3oh = all_lines['CH3OH_422_312']
    nu_h2co = all_lines['H2CO_322_221']
    v_ch3oh = ((nu_ch3oh - nu_h2co)/nu_h2co * constants.c).to(u.km/u.s).value

    hdu = fits.open('APEX_H2CO_322_221_{0}.fits'.format(smooth))[0]
    dv = hdu.header['CDELT3']
    shift = v_ch3oh / dv
    print "dv,shift: ",dv,shift

    mask = fits.getdata('APEX_H2CO_303_202_{0}_mask.fits'.format(smooth)).astype('bool')
    print mask.shape
    newmask = mask*False
    print newmask.shape
    newmask[np.abs(shift):,:,:] = mask[:-np.abs(shift),:,:]
    print newmask.sum()
    hdu.data[newmask] = np.nan
    hdu.writeto('APEX_H2CO_322_221_{0}_CH3OHchomped.fits'.format(smooth), clobber=True)

    hdu.data[True-mask] = np.nan
    hdu.writeto('APEX_H2CO_322_221_{0}_CH3OHchomped_masked.fits'.format(smooth), clobber=True)

    integrate_mask('APEX_H2CO_322_221_{0}_CH3OHchomped'.format(smooth),
                   mask='APEX_H2CO_303_202_{0}_mask.fits'.format(smooth))

def do_mask_ch3oh():
    # spatial smoothing = 2pix
    mask_out_ch3oh('smooth')
    # spatial smoothing = 4pix
    mask_out_ch3oh('vsmooth')

def do_2014(datasets=datasets_2014):
    #datasets = ['E-093.C-0144A.2014APR02/E-093.C-0144A-2014-2014-04-01',
    #            'E-093.C-0144A.2014APR03/E-093.C-0144A-2014-2014-04-02']
    #build_cube_2014('MAP_001', datasets=datasets, scanblsub=True, lowhigh='low')
    #build_cube_2014('MAP_001', datasets=datasets, scanblsub=True, lowhigh='high')
    #build_cube_2014('MAP_001', datasets=datasets, scanblsub=False, lowhigh='high_nosub')

    for dataset in datasets:
        for source in datasets[dataset]:
            build_cube_2014(source, datasets=[dataset], scanblsub=True, lowhigh='low')
            build_cube_2014(source, datasets=[dataset], scanblsub=True, lowhigh='high')


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


def build_cube_2014(sourcename,
                    mergefile=None,
                    lowhigh='low',
                    downsample_factor=8,
                    datapath='/Users/adam/work/h2co/apex/april2014/',
                    outpath='/Users/adam/work/h2co/apex/april2014/',
                    datasets=None,
                    scanblsub=False,
                    verbose=True,
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

    # need two loops to avoid loading too much stuff into memory
    for dataset in datasets:

        apex_filename=datapath+dataset+".apex"
        
        log.info("".join(("Loading data for dataset ",dataset," to filename ",apex_filename,"  t=",str(time.time()-t0))))

        spectra,headers,indices = load_apex_cube(apex_filename, skip_data=False,
                                                 downsample_factor=downsample_factor,
                                                 )

        #if dataset == 'M-091.F-0019-2013-2013-06-13':
        #    tsysrange=[100,260]
        #else:
        #    tsysrange=[100,325]
        tsysrange=[100,325]

        log.info("".join(("Selecting data for dataset ",dataset," to filename ",apex_filename,"  t=",str(time.time()-t0))))

        data, hdrs, gal = select_apex_data(spectra, headers, indices,
                                           sourcename=sourcename,
                                           # NOT ignored, even though it's not used above...
                                           # this is only OK because the bad shapes are from
                                           # Saturn
                                           #shapeselect=4096,
                                           tsysrange=tsysrange,
                                           xtel=xtel,
                                           rchanrange=None,
                                           skip_data=False)

        log.info("".join(("Processing data for dataset ",dataset," to filename ",apex_filename,"  t=",str(time.time()-t0))))

        data, gal, hdrs = process_data(data, gal, hdrs, os.path.join(outpath,
                                                                     dataset)+"_"+xtel,
                                       scanblsub=scanblsub, verbose=verbose,
                                       **kwargs)

        log.info("".join(("Adding data for dataset ",dataset," to filename ",apex_filename,"  t=",str(time.time()-t0))))

        add_apex_data(data, hdrs, gal, cubefilename, retfreq=True,
                      varweight=True,
                      # downsample factor for freqarr
                      )
        # FORCE cleanup
        log.info("".join(("Clearing data for dataset ",dataset," to filename ",apex_filename,"  t=",str(time.time()-t0))))
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
                             nsigma_ignore=1.0, return_mask=False):
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
        The number of standard deviations above which the mean spectrum ought
        to be before it is ignored
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

    if return_mask and automask > 0:
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
                log.info("Masked {0} pixels in scan {1}-{2} ({3}%)".format(nflag,
                                                                           ii, jj,
                                                                           nflag/float(mask_pixels.size),
                                                                           )
                          )

        if mask_pixels is None:
            y = data[ii:jj,:]
        else:
            # mask_pixels is an include mask
            inds = np.arange(data.shape[freqaxis])[mask_pixels]
            y = data[ii:jj,mask_pixels]
            if return_mask:
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

    if return_mask:
        return dsub, np.array(masklist)

    return dsub

def efuncs(arr, return_others=False):
    """
    Determine eigenfunctions of an array for use with
    PCA cleaning
    """
    if hasattr(arr,'filled'):
        arr = arr.filled(0)
    covmat = np.dot(arr.T,arr)
    evals,evects = np.linalg.eig(covmat)
    efuncarr = np.dot(arr,evects)
    if return_others:
        return efuncarr,covmat,evals,evects
    else:
        return efuncarr

def PCA_clean(data,
              smoothing_scale=200.,
              timeaxis=0,
              freqaxis=1,
              ncomponents=3,
             ):
    """
    Remove N PCA components in the time direction

    TODO: speed up by downsampling in TIME as well; we don't expect large
    second-to-second variations
    """

    if freqaxis == 0 and timeaxis == 1:
        data = data.swapaxes(0,1)
    elif freqaxis != 1 or timeaxis != 0:
        raise ValueError("Invalid axis specification.")

    if np.any(np.isnan(data)):
        warnings.warn("There were NaNs in the PCA-target data")
        data = np.nan_to_num(data)

    log.info(("PCA cleaning an image with size {0},"
              " which will downsample to {1}").format(data.shape,
                                                      (data.shape[0],
                                                       data.shape[1]/(smoothing_scale/5))))

    sm_data = filters.gaussian_filter1d(data, smoothing_scale,
                                        axis=1, mode='mirror')

    efuncarr,covmat,evals,evects = efuncs(sm_data[:,::smoothing_scale/5].T,
                                          return_others=True)

    # Zero-out the components we want to keep
    efuncarr[:,ncomponents:] = 0

    to_subtract = np.inner(efuncarr,evects).T

    ifunc = interpolate.interp1d(np.arange(to_subtract.shape[1]),
                                 to_subtract,
                                 axis=1)
    to_subtract = ifunc(np.linspace(0, to_subtract.shape[1]-1, data.shape[1]))

    dsub = data - to_subtract

    if freqaxis == 0 and timeaxis == 1:
        dsub = dsub.swapaxes(0,1)

    return dsub


def demo_parameters_of_blsubbing(data, scans, mask):
    """
    These figures are used to justify the choice of parameters in
    subtract_scan_lienar_fit:
        automask=2, smooth_all=True

    spectra,headers,indices,data,hdrs,gal,dspecsub,dmeansub,freq,mask = debug_and_load()
    """

    data_diagplot(data, 'blsub_demo_raw')

    dsub = subtract_scan_linear_fit(data, scans, verbose=True, automask=True,
                                    smooth_all=True)
    data_diagplot(dsub, 'blsub_demo_automask_smoothall')

    dsub = subtract_scan_linear_fit(data, scans, verbose=True, automask=True,
                                    smooth_all=False)
    data_diagplot(dsub, 'blsub_demo_automask')

    dsub = subtract_scan_linear_fit(data, scans, verbose=True, automask=2,
                                    smooth_all=False)
    data_diagplot(dsub, 'blsub_demo_automask_smooth2')

    dsub = subtract_scan_linear_fit(data, scans, mask_pixels=mask,
                                    verbose=True, automask=False,
                                    smooth_all=False)
    data_diagplot(dsub, 'blsub_demo_linemask')

    dsub = subtract_scan_linear_fit(data, scans, verbose=True, automask=2,
                                    smooth_all=True)
    data_diagplot(dsub, 'blsub_demo_automask_smooth2_smoothall')

#make_high_mergecube()
#make_low_mergecube()
