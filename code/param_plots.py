import paths
import os
import pylab as pl
import numpy as np
from astropy import units as u
from pyspeckit_fitting import (texgrid303, taugrid303, texgrid321, taugrid321,
                               texgrid322, taugrid322, hdr)

from h2co_modeling import grid_fitter
from astropy import table
from scipy.ndimage.interpolation import map_coordinates

pl.rcParams['font.size'] = 16.0

Tbackground=2.73
tline303a = (1.0-np.exp(-np.array(taugrid303)))*(texgrid303-Tbackground)
tline321a = (1.0-np.exp(-np.array(taugrid321)))*(texgrid321-Tbackground)
tline322a = (1.0-np.exp(-np.array(taugrid322)))*(texgrid322-Tbackground)

zinds,yinds,xinds = np.indices(tline303a.shape)
upsample_factor = np.array([5,12.5,12.5], dtype='float')
uzinds,uyinds,uxinds = upsinds = np.indices([x*us
                                             for x,us in zip(tline303a.shape,
                                                             upsample_factor)],
                                           dtype='float')
tline303 = map_coordinates(tline303a,
                           upsinds/upsample_factor[:,None,None,None],
                           mode='nearest')
tline321 = map_coordinates(tline321a,
                           upsinds/upsample_factor[:,None,None,None],
                           mode='nearest')
densityarr = ((uxinds + hdr['CRPIX1']-1)*hdr['CDELT1'] /
              float(upsample_factor[2])+hdr['CRVAL1']) # log density
columnarr  = ((uyinds + hdr['CRPIX2']-1)*hdr['CDELT2'] /
              float(upsample_factor[1])+hdr['CRVAL2']) # log column
temparr    = ((uzinds + hdr['CRPIX3']-1)*hdr['CDELT3'] /
              float(upsample_factor[0])+hdr['CRVAL3']) # lin temperature
drange = [densityarr.min(), densityarr.max()]
crange = [columnarr.min(), columnarr.max()]
trange = [temparr.min(), temparr.max()]
darr = densityarr[0,0,:]
carr = columnarr[0,:,0]
tarr = temparr[:,0,0]

# While the individual lines are subject to filling factor uncertainties, the
# ratio is not.
modelratio = tline321/tline303

fittable = table.Table.read("fitted_line_parameters.ipac", format='ascii.ipac')

if not os.path.exists(paths.fpath('param_fits')):
    os.makedirs(paths.fpath('param_fits'))

arbitrary_scale = 1
nlevs=5

for row in fittable:
    logh2column = np.log10(row['higalcolumndens'])
    elogh2column = 1.0
    linewidth = row['width_0']
    elinewidth = row['ewidth_0']

    par1 = row['ampH2CO_0']
    epar1 = row['eampH2CO_0']*arbitrary_scale
    par2 = row['ampH2CO_0']*row['h2coratio_0']
    epar2 = row['ampH2CO_0']*row['eh2coratio_0']*arbitrary_scale
    #match,indbest,chi2b = grid_fitter.grid_2p_getmatch(par1, epar1, tline303,
    #                                                   par2, epar2, tline321)
    ratio = row['h2coratio_0']
    eratio = row['eh2coratio_0']
    match,indbest,chi2a = grid_fitter.grid_getmatch(ratio, eratio, modelratio)

    # We can impose a "loose" abundance constraint
    # Given that we know the H2 density, and the line width is ~5-10 km/s,
    # abundance = column / pc / density
    # We'll say abundance = 1.2e9 with error 0.6e9
    # Or, log(abundance) = log(1.2e9) +/- 1
    logabundance = np.log10(1.2e-9)
    elogabundance = 1.0
    model_logabundance = np.log10(10**columnarr / u.pc.to(u.cm) / 10**densityarr)
    chi2X = ((model_logabundance-logabundance)/elogabundance)**2

    # Combined abundance + total column constraint
    # N(H2CO) * dv * X = N(H2)
    # We are effectively ignoring errors in the linewidth here:
    h2fromh2co = np.log10(10**columnarr * (np.sqrt(np.pi) * linewidth) * 10**logabundance)
    chi2_h2 = ((h2fromh2co-logh2column)/elogh2column)**2

    # Even though the lines are subject to filling-factor uncertainty, we can
    # set a *minimum* brightness in the models.  Given that we observe a line
    # brightness T_A, the true brightness is T_B = T_A/ff, where ff<1 by
    # definition
    # We therefore *increase* the chi^2 value wherever the model is fainter
    # than the line, enforcing a soft lower limit
    chi2_1 = ((tline303 - par1)/epar1)**2 * (tline303 < par1)
    chi2_2 = ((tline321 - par2)/epar2)**2 * (tline321 < par2)
    chi2b = chi2a + chi2_1 + chi2_2 + chi2X + chi2_h2
    match = chi2b < 1
    indbest,match = grid_fitter.getmatch(chi2b, match)

    sh = match.shape
    (zz,yy,xx) = np.unravel_index(indbest, sh)

    pl.figure(1)
    pl.clf()
    pl.subplot(2,3,1)
    pl.imshow(tline303[zz,:,:], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=drange+crange)
    pl.contour(darr, carr, chi2b[zz,:,:], levels=chi2b.min()+np.arange(nlevs))
    pl.ylabel('Column')
    pl.xlabel('Density')

    pl.subplot(2,3,2)
    pl.imshow(tline303[:,yy,:], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(drange)/np.diff(trange),
              extent=drange+trange)
    pl.contour(darr, tarr, chi2b[:,yy,:], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel('Density')
    pl.ylabel('Temperature')

    pl.subplot(2,3,3)
    pl.imshow(tline303[:,:,xx], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(crange)/np.diff(trange),
              extent=crange+trange)
    pl.contour(carr, tarr, chi2b[:,:,xx], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel('Column')
    pl.ylabel('Temperature')

    pl.subplot(2,3,4)
    pl.imshow(tline321[zz,:,:], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=drange+crange)
    pl.contour(darr, carr, chi2b[zz,:,:], levels=chi2b.min()+np.arange(nlevs))
    pl.ylabel('Column')
    pl.xlabel('Density')

    pl.subplot(2,3,5)
    pl.imshow(tline321[:,yy,:], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(drange)/np.diff(trange),
              extent=drange+trange)
    pl.contour(darr, tarr, chi2b[:,yy,:], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel('Density')
    pl.ylabel('Temperature')

    pl.subplot(2,3,6)
    pl.imshow(tline321[:,:,xx], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(crange)/np.diff(trange),
              extent=crange+trange)
    pl.contour(carr, tarr,  chi2b[:,:,xx], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel('Column')
    pl.ylabel('Temperature')

    pl.suptitle(row['Source_Name'])

    pl.savefig(paths.fpath('param_fits/{name}_h2coratio.pdf'.format(name=row['Source_Name'])))

    pl.figure(2)
    pl.clf()
    pl.subplot(2,3,1)
    yi, xi = np.indices(tline303.shape[1:])
    inds = [chi2b.argmin(axis=0), yi, xi]
    pl.imshow(tline303[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=drange+crange)
    pl.contour(darr, carr, chi2b.min(axis=0), levels=chi2b.min()+np.arange(nlevs))
    pl.ylabel('Column')
    pl.xlabel('Density')

    pl.subplot(2,3,2)
    zi, xi = np.indices([tline303.shape[0], tline303.shape[2],])
    inds = [zi, chi2b.argmin(axis=1), xi]
    pl.imshow(tline303[:,yy,:], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(drange)/np.diff(trange),
              extent=drange+trange)
    pl.contour(darr, tarr, chi2b.min(axis=1), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel('Density')
    pl.ylabel('Temperature')

    pl.subplot(2,3,3)
    zi, yi = np.indices([tline303.shape[0], tline303.shape[2],])
    inds = [zi, yi, chi2b.argmin(axis=2)]
    pl.imshow(tline303[:,:,xx], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(crange)/np.diff(trange),
              extent=crange+trange)
    pl.contour(carr, tarr, chi2b.min(axis=2), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel('Column')
    pl.ylabel('Temperature')

    pl.subplot(2,3,4)
    yi, xi = np.indices(tline303.shape[1:])
    inds = [chi2b.argmin(axis=0), yi, xi]
    pl.imshow(tline321[inds], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=drange+crange)
    pl.contour(darr, carr, chi2b.min(axis=0), levels=chi2b.min()+np.arange(nlevs))
    pl.ylabel('Column')
    pl.xlabel('Density')

    pl.subplot(2,3,5)
    zi, xi = np.indices([tline303.shape[0], tline303.shape[2],])
    inds = [zi, chi2b.argmin(axis=1), xi]
    pl.imshow(tline321[inds], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(drange)/np.diff(trange),
              extent=drange+trange)
    pl.contour(darr, tarr, chi2b.min(axis=1), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel('Density')
    pl.ylabel('Temperature')

    pl.subplot(2,3,6)
    zi, yi = np.indices([tline303.shape[0], tline303.shape[2],])
    inds = [zi, yi, chi2b.argmin(axis=2)]
    pl.imshow(tline321[inds], cmap=pl.cm.bone, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(crange)/np.diff(trange),
              extent=crange+trange)
    pl.contour(carr, tarr,  chi2b.min(axis=2), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel('Column')
    pl.ylabel('Temperature')

    pl.suptitle(row['Source_Name'])

    pl.savefig(paths.fpath('param_fits/{name}_h2coratio_minaxis.pdf'.format(name=row['Source_Name'])))
    #break
pl.show()
