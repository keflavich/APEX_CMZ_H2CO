import paths
import os
import pylab as pl
import numpy as np
from astropy import log
from astropy import units as u
from paths import analysispath
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

fittable = table.Table.read(os.paths.join(analysispath,
                                          "fitted_line_parameters.ipac"),
                            format='ascii.ipac')
fittable.add_columns([table.Column(name=name, dtype='float', length=len(fittable))
                      for name in ['temperature_chi2','tmin1sig_chi2','tmax1sig_chi2',
                                   'column_chi2','cmin1sig_chi2','cmax1sig_chi2',
                                   'density_chi2','dmin1sig_chi2','dmax1sig_chi2',]])

if not os.path.exists(paths.fpath('param_fits')):
    os.makedirs(paths.fpath('param_fits'))

arbitrary_scale = 1
nlevs=5

density_label = 'Density $n(\mathrm{H}_2)$ [log cm$^{-3}$]'
column_label = 'p-H$_2$CO [log cm$^{-2}$/(km s$^{-1}$ pc)]'
temperature_label = 'Temperature (K)'

for row in fittable:
    log.info(row['Source_Name'])
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
    match,indbest,chi2r = grid_fitter.grid_getmatch(ratio, eratio, modelratio)

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
    h2fromh2co = np.log10(10**columnarr * (np.sqrt(np.pi) * linewidth) / 10**logabundance)
    chi2_h2 = ((h2fromh2co-logh2column)/elogh2column)**2

    # Even though the lines are subject to filling-factor uncertainty, we can
    # set a *minimum* brightness in the models.  Given that we observe a line
    # brightness T_A, the true brightness is T_B = T_A/ff, where ff<1 by
    # definition
    # We therefore *increase* the chi^2 value wherever the model is fainter
    # than the line, enforcing a soft lower limit
    chi2_1 = ((tline303 - par1)/epar1)**2 * (tline303 < par1)
    chi2_2 = ((tline321 - par2)/epar2)**2 * (tline321 < par2)
    chi2_ff = chi2_1+chi2_2
    chi2b = chi2r + chi2_ff + chi2X + chi2_h2
    match = chi2b < 1
    indbest,match = grid_fitter.getmatch(chi2b, match)

    sh = match.shape
    (zz,yy,xx) = np.unravel_index(indbest, sh)

    fig1 = pl.figure(1)
    fig1.clf()

    vmin = np.max([tline303.min(), 0.1])
    vmax = np.min([tline303.max(), par1+10])
    ax1 = pl.subplot(2,3,1)
    im1 = pl.imshow(tline303[zz,:,:], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=drange+crange, vmin=vmin, vmax=vmax)
    pl.contour(darr, carr, chi2b[zz,:,:], levels=chi2b.min()+np.arange(nlevs))
    pl.ylabel(column_label)
    pl.xlabel(density_label)

    ax2 = pl.subplot(2,3,2)
    im2 = pl.imshow(tline303[:,yy,:], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(drange)/np.diff(trange),
              extent=drange+trange, vmin=vmin, vmax=vmax)
    pl.contour(darr, tarr, chi2b[:,yy,:], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(density_label)
    pl.ylabel(temperature_label)
    #ax2.set_title("p-H$_2$CO $3_{0,3}-2_{0,2}$")

    ax3 = pl.subplot(2,3,3)
    im3 = pl.imshow(tline303[:,:,xx], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(crange)/np.diff(trange),
              extent=crange+trange, vmin=vmin, vmax=vmax)
    pl.contour(carr, tarr, chi2b[:,:,xx], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label)
    ax3.xaxis.set_ticks(np.arange(carr.min(), carr.max()))
    pl.ylabel(temperature_label)
    cax = fig1.add_axes([0.91,0.55,0.02,0.35])
    cb = fig1.colorbar(mappable=im3, cax=cax, ax=ax2)
    cb.set_label("$T_B$ (p-H$_2$CO $3_{0,3}-2_{0,2}$)")

    vmin = np.max([tline321.min(), 0.1])
    vmax = np.min([tline321.max(), par2+10])
    ax4 = pl.subplot(2,3,4)
    pl.imshow(tline321[zz,:,:], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=drange+crange, vmin=vmin, vmax=vmax)
    pl.contour(darr, carr, chi2b[zz,:,:], levels=chi2b.min()+np.arange(nlevs))
    pl.ylabel(column_label)
    pl.xlabel(density_label)

    ax5 = pl.subplot(2,3,5)
    im5 = pl.imshow(tline321[:,yy,:], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(drange)/np.diff(trange),
              extent=drange+trange, vmin=vmin, vmax=vmax)
    pl.contour(darr, tarr, chi2b[:,yy,:], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(density_label)
    pl.ylabel(temperature_label)
    #ax5.set_title("p-H$_2$CO $3_{2,1}-2_{2,0}$")

    ax6 = pl.subplot(2,3,6)
    im6 = pl.imshow(tline321[:,:,xx], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(crange)/np.diff(trange),
              extent=crange+trange, vmin=vmin, vmax=vmax)
    pl.contour(carr, tarr,  chi2b[:,:,xx], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label)
    ax6.xaxis.set_ticks(np.arange(carr.min(), carr.max()))
    pl.ylabel(temperature_label)
    cax = fig1.add_axes([0.91,0.1,0.02,0.35])
    cb = fig1.colorbar(mappable=im6, cax=cax, ax=ax5)
    cb.set_label("$T_B$ (p-H$_2$CO $3_{2,1}-2_{2,0}$)")

    pl.suptitle(row['Source_Name'])
    pl.subplots_adjust(wspace=0.33, hspace=0.22, left=0.1)

    pl.savefig(paths.fpath('param_fits/{name}_h2coratio.pdf'.format(name=row['Source_Name'])), bbox_inches='tight')

    fig2 = pl.figure(2)
    fig2.clf()
    ax1 = pl.subplot(2,3,1)
    yi, xi = np.indices(tline303.shape[1:])
    inds = [chi2b.argmin(axis=0), yi, xi]
    # The background from taking the min-chi^2 along each axis is too ugly and
    # hard to explain: revert to using a *slice* for a background but a chi^2
    # *projection* for the contours
    inds = [zz, slice(None), slice(None)]
    pl.imshow(tline303[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=drange+crange, vmin=vmin, vmax=vmax)
    pl.contour(darr, carr, chi2b.min(axis=0), levels=chi2b.min()+np.arange(nlevs))
    pl.ylabel(column_label)
    pl.xlabel(density_label)

    ax2 = pl.subplot(2,3,2)
    zi, xi = np.indices([tline303.shape[0], tline303.shape[2],])
    inds = [zi, chi2b.argmin(axis=1), xi]
    inds = [slice(None), yy, slice(None)]
    pl.imshow(tline303[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(drange)/np.diff(trange),
              extent=drange+trange, vmin=vmin, vmax=vmax)
    pl.contour(darr, tarr, chi2b.min(axis=1), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(density_label)
    pl.ylabel(temperature_label)

    ax3 = pl.subplot(2,3,3)
    zi, yi = np.indices([tline303.shape[0], tline303.shape[2],])
    inds = [zi, yi, chi2b.argmin(axis=2)]
    inds = [slice(None), slice(None), xx]
    pl.imshow(tline303[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(crange)/np.diff(trange),
              extent=crange+trange, vmin=vmin, vmax=vmax)
    pl.contour(carr, tarr, chi2b.min(axis=2), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label)
    ax3.xaxis.set_ticks(np.arange(carr.min(), carr.max()))
    pl.ylabel(temperature_label)
    cax = fig2.add_axes([0.91,0.55,0.02,0.35])
    cb = fig2.colorbar(mappable=im3, cax=cax, ax=ax2)
    cb.set_label("$T_B$ (p-H$_2$CO $3_{0,3}-2_{0,2}$)")

    ax4 = pl.subplot(2,3,4)
    yi, xi = np.indices(tline303.shape[1:])
    inds = [chi2b.argmin(axis=0), yi, xi]
    inds = [zz, slice(None), slice(None)]
    pl.imshow(tline321[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=drange+crange, vmin=vmin, vmax=vmax)
    pl.contour(darr, carr, chi2b.min(axis=0), levels=chi2b.min()+np.arange(nlevs))
    pl.ylabel(column_label)
    pl.xlabel(density_label)

    ax5 = pl.subplot(2,3,5)
    zi, xi = np.indices([tline303.shape[0], tline303.shape[2],])
    inds = [zi, chi2b.argmin(axis=1), xi]
    inds = [slice(None), yy, slice(None)]
    pl.imshow(tline321[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(drange)/np.diff(trange),
              extent=drange+trange, vmin=vmin, vmax=vmax)
    pl.contour(darr, tarr, chi2b.min(axis=1), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(density_label)
    pl.ylabel(temperature_label)

    ax6 = pl.subplot(2,3,6)
    zi, yi = np.indices([tline303.shape[0], tline303.shape[2],])
    inds = [zi, yi, chi2b.argmin(axis=2)]
    inds = [slice(None), slice(None), xx]
    im6 = pl.imshow(tline321[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(crange)/np.diff(trange),
              extent=crange+trange, vmin=vmin, vmax=vmax)
    pl.contour(carr, tarr,  chi2b.min(axis=2), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label)
    ax6.xaxis.set_ticks(np.arange(carr.min(), carr.max()))
    pl.ylabel(temperature_label)
    cax = fig2.add_axes([0.91,0.1,0.02,0.35])
    cb = fig2.colorbar(mappable=im6, cax=cax, ax=ax5)
    cb.set_label("$T_B$ (p-H$_2$CO $3_{2,1}-2_{2,0}$)")

    pl.suptitle(row['Source_Name'])
    pl.subplots_adjust(wspace=0.33, left=0.1, hspace=0.22)

    pl.savefig(paths.fpath('param_fits/{name}_h2coratio_minaxis.pdf'.format(name=row['Source_Name'])), bbox_inches='tight')

    # Show the constraints provided by individual parameters
    pl.figure(3)
    pl.clf()
    # chi2b = chi2r + chi2_1 + chi2_2 + chi2X + chi2_h2
    ax1 = pl.subplot(2,3,1)
    pl.contourf(darr, tarr, chi2r.min(axis=1), levels=chi2r.min()+np.arange(nlevs), alpha=0.5)
    pl.xlabel(density_label)
    pl.ylabel(temperature_label)
    pl.title("Ratio $3_{0,3}-2_{0,2}/3_{2,1}-2_{2,0}$")
    ax2 = pl.subplot(2,3,2)
    pl.contourf(darr, tarr, chi2r.min(axis=0), levels=chi2r.min()+np.arange(nlevs), alpha=0.5)
    pl.xlabel(density_label)
    pl.ylabel(column_label)
    pl.title("Ratio $3_{0,3}-2_{0,2}/3_{2,1}-2_{2,0}$")
    ax4 = pl.subplot(2,3,6)
    pl.contourf(darr, carr, chi2X.min(axis=0), levels=chi2r.min()+np.arange(nlevs), alpha=0.5)
    pl.ylabel(column_label)
    pl.xlabel(density_label)
    pl.title("log(p-H$_2$CO/H$_2$) $= {0:0.1f}\pm{1:0.1f}$".format(logabundance, elogabundance))
    ax3 = pl.subplot(2,3,3)
    pl.contourf(darr, carr, chi2_h2.min(axis=0), levels=chi2_h2.min()+np.arange(nlevs), alpha=0.5)
    pl.xlabel(density_label)
    pl.ylabel(column_label)
    pl.title("Total log$(N(\\mathrm{{H}}_2)) = {0:0.1f}\pm{1:0.1f}$".format(logh2column,
                                                                                    elogh2column))
    ax5 = pl.subplot(2,3,5)
    pl.contourf(darr, carr, (chi2_ff.min(axis=0)),
                levels=chi2_ff.min()+np.arange(nlevs), alpha=0.5)
    pl.xlabel(density_label)
    pl.ylabel(column_label)
    pl.title("Line Brightness + $ff\leq1$")
    ax6 = pl.subplot(2,3,4)
    pl.contourf(darr, tarr, (chi2_ff.min(axis=1)),
                levels=chi2_ff.min()+np.arange(nlevs), alpha=0.5)
    pl.xlabel(density_label)
    pl.ylabel(temperature_label)
    pl.title("Line Brightness + $ff\leq1$")

    pl.subplots_adjust(wspace=0.4, hspace=0.4)
    pl.savefig(paths.fpath('param_fits/{name}_parameter_constraints.pdf'.format(name=row['Source_Name'])), bbox_inches='tight')
    #break

    deltachi2b = (chi2b-chi2b.min())
    for parname,pararr in zip(('temperature','column','density'),
                              (temparr,columnarr,densityarr)):
        row['{0}_chi2'.format(parname)] = pararr.flat[indbest]
        OK = deltachi2b<1
        if np.count_nonzero(OK) > 0:
            row['{0:1.1s}min1sig_chi2'.format(parname)] = pararr[OK].min()
            row['{0:1.1s}max1sig_chi2'.format(parname)] = pararr[OK].max()
        else:
            row['{0:1.1s}min1sig_chi2'.format(parname)] = np.nan
            row['{0:1.1s}max1sig_chi2'.format(parname)] = np.nan


fittable.write(os.path.join(analysispath,
                            'fitted_line_parameters_Chi2Constraints.ipac'),
               format='ascii.ipac')

pl.show()
