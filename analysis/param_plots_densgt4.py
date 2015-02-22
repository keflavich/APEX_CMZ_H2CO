"""
ADDED CONSTRAINT: density > 10^4

For each of the fitted spectra from individual_spectra.py, use the fitted ratio
(and appropriate ancillary data, e.g. h2 column) to derive the best fit
temperature, etc.

Store them in an output table and plot by calling 
execfile(paths.pcpath('parameter_comparisons.py'))
"""
import paths
import os
import pylab as pl
import numpy as np
from astropy import log
from astropy import units as u
from astropy import constants
from astropy import coordinates
from paths import analysispath
from pyspeckit_fitting import (texgrid303, taugrid303, texgrid321, taugrid321,
                               texgrid322, taugrid322, hdr)

from constrain_parameters import paraH2COmodel
from h2co_modeling import grid_fitter
from astropy import table
import heating


pl.rcParams['font.size'] = 16.0

# mf means modelfitter
mf = paraH2COmodel()

fittable = table.Table.read(os.path.join(analysispath,
                                          "fitted_line_parameters.ipac"),
                            format='ascii.ipac')
fittable.add_columns([table.Column(name=name, dtype='float', length=len(fittable))
                      for name in ['temperature_chi2','tmin1sig_chi2','tmax1sig_chi2',
                                   'column_chi2','cmin1sig_chi2','cmax1sig_chi2',
                                   'density_chi2','dmin1sig_chi2','dmax1sig_chi2',
                                   'logh2column','elogh2column',
                                   'logabundance','elogabundance',
                                   'tkin_turb', 'reff_pc',
                                  ]])

if not os.path.exists(paths.fpath('param_fits')):
    os.makedirs(paths.fpath('param_fits'))

nlevs=5

density_label = 'Density $n(\mathrm{H}_2)$ [log cm$^{-3}$]'
column_label = 'Column p-H$_2$CO [log cm$^{-2}$/(km s$^{-1}$ pc)]'
density_label_short = "$n(\mathrm{H}_2) (\mathrm{cm}^{-3})$"
column_label_short = "$N(\mathrm{p-H}_2\mathrm{CO}) (\mathrm{cm}^{-2})$"
temperature_label = 'Temperature (K)'
prevname = ''
num = 0

for row in fittable:
    if row['Source_Name'] == prevname:
        num += 1
    else:
        num = 0
        prevname = row['Source_Name']
    log.info("Fitting {0}_{1}".format(row['Source_Name'],num))
    logh2column = np.log10(row['higalcolumndens'])
    elogh2column = 1.0
    linewidth = row['spline_width']
    elinewidth = row['espline_width']

    par1 = row['spline_ampH2CO']
    epar1 = row['espline_ampH2CO']
    par2 = row['spline_ampH2CO']*row['spline_h2coratio321303']
    epar2 = row['spline_ampH2CO']*row['espline_h2coratio321303']
    #match,indbest,chi2b = grid_fitter.grid_2p_getmatch(par1, epar1, tline303,
    #                                                   par2, epar2, tline321)
    ratio = row['spline_h2coratio321303']
    eratio = row['espline_h2coratio321303']
    ratio2 = row['spline_h2coratio322321']
    eratio2 = row['espline_h2coratio322321']

    # We can impose a "loose" abundance constraint
    # Given that we know the H2 density, and the line width is ~5-10 km/s,
    # abundance = column / pc / density
    # We'll say abundance = 1.2e9 with error 0.6e9
    # Or, log(abundance) = log(1.2e9) +/- 1
    logabundance = np.log10(1.2e-9)
    elogabundance = 1.0

    # Constraint from density
    r_deg = (row['area']/np.pi)**0.5 * u.deg
    reff = (r_deg*(8.5*u.kpc)).to(u.pc, u.dimensionless_angles())
    mass = ((10**logh2column*u.cm**-2)*np.pi*reff**2*2.8*constants.m_p).to(u.M_sun)
    density = (mass/(4/3.*np.pi*reff**3)/constants.m_p/2.8).to(u.cm**-3)
    mindens = np.log10(density.value)
    if mindens < 3:
        mindens = 3

    # Combined abundance + total column constraint
    # N(H2CO) * dv * X = N(H2)
    # We are effectively ignoring errors in the linewidth here:
    # (noop - see chi2_h2)

    # Even though the lines are subject to filling-factor uncertainty, we can
    # set a *minimum* brightness in the models.  Given that we observe a line
    # brightness T_A, the true brightness is T_B = T_A/ff, where ff<1 by
    # definition
    # We therefore *increase* the chi^2 value wherever the model is fainter
    # than the line, enforcing a soft lower limit

    mf.set_constraints(ratio303321=ratio, eratio303321=eratio,
                       ratio321322=ratio2, eratio321322=eratio2,
                       logh2column=logh2column, elogh2column=elogh2column,
                       logabundance=logabundance, elogabundance=elogabundance,
                       taline303=par1, etaline303=epar1,
                       taline321=par2, etaline321=epar2,
                       mindens=mindens,
                       linewidth=linewidth)


    chi2r = mf.chi2_r303321
    chi2r2 = mf.chi2_r321322
    chi2_h2 = mf.chi2_h2
    chi2X = mf.chi2_X
    chi2_1 = mf.chi2_ff1
    chi2_2 = mf.chi2_ff2
    chi2_ff = chi2_1+chi2_2
    chi2_dens = mf.chi2_dens
    chi2b = chi2r + chi2_ff + chi2X + chi2_h2 + chi2_dens

    match = chi2b < 1
    indbest,match = grid_fitter.getmatch(chi2b, match)

    sh = match.shape
    (zz,yy,xx) = np.unravel_index(indbest, sh)

    fig1 = pl.figure(1, figsize=(12,12))
    fig1.clf()

    vmin = np.max([mf.tline303.min(), 0.1])
    vmax = np.min([mf.tline303.max(), par1+10])
    ax1 = pl.subplot(3,3,1)
    im1 = pl.imshow(mf.tline303[zz,:,:], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=mf.crange+mf.drange, vmin=vmin, vmax=vmax)
    pl.contour(mf.carr, mf.darr, chi2b[zz,:,:], levels=chi2b.min()+np.arange(nlevs))
    # only label middle row
    #pl.ylabel(column_label)
    pl.xlabel(column_label_short)
    ax1.set_xticks(np.arange(11,16))

    ax2 = pl.subplot(3,3,2)
    im2 = pl.imshow(mf.tline303[:,yy,:], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(mf.crange)/np.diff(mf.trange),
              extent=mf.crange+mf.trange, vmin=vmin, vmax=vmax)
    pl.contour(mf.carr, mf.tarr, chi2b[:,yy,:], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label_short)
    ax2.set_xticks(np.arange(11,16))
    pl.ylabel(temperature_label)
    #ax2.set_title("p-H$_2$CO $3_{0,3}-2_{0,2}$")

    ax3 = pl.subplot(3,3,3)
    im3 = pl.imshow(mf.tline303[:,:,xx], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(mf.drange)/np.diff(mf.trange),
              extent=mf.drange+mf.trange, vmin=vmin, vmax=vmax)
    pl.contour(mf.darr, mf.tarr, chi2b[:,:,xx], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(density_label_short)
    ax3.xaxis.set_ticks(np.arange(mf.darr.min(), mf.darr.max()))
    pl.ylabel(temperature_label)
    cax = fig1.add_axes([0.91,0.68,0.02,0.22])
    cb = fig1.colorbar(mappable=im3, cax=cax, ax=ax2)
    cb.set_label("$T_B$ (p-H$_2$CO $3_{0,3}-2_{0,2}$)")

    vmin = np.max([mf.tline321.min(), 0.1])
    vmax = np.min([mf.tline321.max(), par2+10])
    ax4 = pl.subplot(3,3,4)
    pl.imshow(mf.tline321[zz,:,:], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=mf.crange+mf.drange, vmin=vmin, vmax=vmax)
    pl.contour(mf.carr, mf.darr, chi2b[zz,:,:], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label_short)
    ax4.set_xticks(np.arange(11,16))
    pl.ylabel(density_label_short)

    ax5 = pl.subplot(3,3,5)
    im5 = pl.imshow(mf.tline321[:,yy,:], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(mf.crange)/np.diff(mf.trange),
              extent=mf.crange+mf.trange, vmin=vmin, vmax=vmax)
    pl.contour(mf.carr, mf.tarr, chi2b[:,yy,:], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label_short)
    ax5.set_xticks(np.arange(11,16))
    pl.ylabel(temperature_label)
    #ax5.set_title("p-H$_2$CO $3_{2,1}-2_{2,0}$")

    ax6 = pl.subplot(3,3,6)
    im6 = pl.imshow(mf.tline321[:,:,xx], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(mf.drange)/np.diff(mf.trange),
              extent=mf.drange+mf.trange, vmin=vmin, vmax=vmax)
    pl.contour(mf.darr, mf.tarr,  chi2b[:,:,xx], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(density_label_short)
    ax6.xaxis.set_ticks(np.arange(mf.darr.min(), mf.darr.max()))
    pl.ylabel(temperature_label)
    cax = fig1.add_axes([0.91,0.40,0.02,0.22])
    cb = fig1.colorbar(mappable=im6, cax=cax, ax=ax5)
    cb.set_label("$T_B$ (p-H$_2$CO $3_{2,1}-2_{2,0}$)")

    vminr = 0.05
    vmaxr = 0.7
    ax7 = pl.subplot(3,3,7)
    im7 = ax7.imshow(mf.modelratio1[zz,:,:], cmap=pl.cm.bone_r,
                     interpolation='spline36',
                     #norm=pl.matplotlib.colors.LogNorm(),
                     extent=mf.crange+mf.drange, vmin=vminr, vmax=vmaxr)
    pl.contour(mf.carr, mf.darr, chi2b[zz,:,:], levels=chi2b.min()+np.arange(nlevs))
    # only label middle row
    #pl.ylabel(column_label_short)
    pl.xlabel(column_label_short)
    ax7.set_xticks(np.arange(11,16))

    ax8 = pl.subplot(3,3,8)
    im8 = ax8.imshow(mf.modelratio1[:,yy,:], cmap=pl.cm.bone_r,
                     interpolation='spline36',
                     #norm=pl.matplotlib.colors.LogNorm(),
                      aspect=np.diff(mf.crange)/np.diff(mf.trange),
                      extent=mf.crange+mf.trange, vmin=vminr, vmax=vmaxr)
    pl.contour(mf.carr, mf.tarr, chi2b[:,yy,:], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label_short)
    ax8.set_xticks(np.arange(11,16))
    pl.ylabel(temperature_label)
    #ax2.set_title("p-H$_2$CO $3_{0,3}-2_{0,2}$")

    ax9 = pl.subplot(3,3,9)
    im9 = ax9.imshow(mf.modelratio1[:,:,xx], cmap=pl.cm.bone_r,
                     interpolation='spline36',
                     #norm=pl.matplotlib.colors.LogNorm(),
                     aspect=np.diff(mf.drange)/np.diff(mf.trange),
                     extent=mf.drange+mf.trange, vmin=vminr, vmax=vmaxr)
    pl.contour(mf.darr, mf.tarr, chi2b[:,:,xx], levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(density_label_short)
    ax9.xaxis.set_ticks(np.arange(mf.darr.min(), mf.darr.max()))
    pl.ylabel(temperature_label)
    cax3 = fig1.add_axes([0.91,0.1,0.02,0.22])
    cb = fig1.colorbar(mappable=im8, cax=cax3, ax=ax8)
    cb.set_label("$3_{2,1}-2_{2,0}$ / $3_{0,3}-2_{0,2}$")


    pl.suptitle(row['Source_Name'])
    pl.subplots_adjust(wspace=0.33, hspace=0.25, left=0.1)

    pl.savefig(paths.fpath('param_fits/{name}_{num}_h2coratio.pdf'.format(name=row['Source_Name'],
                                                                          num=num)), bbox_inches='tight')

    fig2 = pl.figure(2)
    fig2.clf()
    ax1 = pl.subplot(2,3,1)
    yi, xi = np.indices(mf.tline303.shape[1:])
    inds = [chi2b.argmin(axis=0), yi, xi]
    # The background from taking the min-chi^2 along each axis is too ugly and
    # hard to explain: revert to using a *slice* for a background but a chi^2
    # *projection* for the contours
    inds = [zz, slice(None), slice(None)]
    pl.imshow(mf.tline303[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=mf.crange+mf.drange, vmin=vmin, vmax=vmax)
    pl.contour(mf.carr, mf.darr, chi2b.min(axis=0), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label_short)
    ax1.set_xticks(np.arange(11,16))
    pl.ylabel(density_label_short)

    ax2 = pl.subplot(2,3,2)
    zi, xi = np.indices([mf.tline303.shape[0], mf.tline303.shape[2],])
    inds = [zi, chi2b.argmin(axis=1), xi]
    inds = [slice(None), yy, slice(None)]
    pl.imshow(mf.tline303[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(mf.crange)/np.diff(mf.trange),
              extent=mf.crange+mf.trange, vmin=vmin, vmax=vmax)
    pl.contour(mf.carr, mf.tarr, chi2b.min(axis=1), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label_short)
    ax2.set_xticks(np.arange(11,16))
    pl.ylabel(temperature_label)

    ax3 = pl.subplot(2,3,3)
    zi, yi = np.indices([mf.tline303.shape[0], mf.tline303.shape[2],])
    inds = [zi, yi, chi2b.argmin(axis=2)]
    inds = [slice(None), slice(None), xx]
    pl.imshow(mf.tline303[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(mf.drange)/np.diff(mf.trange),
              extent=mf.drange+mf.trange, vmin=vmin, vmax=vmax)
    pl.contour(mf.darr, mf.tarr, chi2b.min(axis=2), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(density_label_short)
    ax3.xaxis.set_ticks(np.arange(mf.darr.min(), mf.darr.max()))
    pl.ylabel(temperature_label)
    cax = fig2.add_axes([0.91,0.55,0.02,0.35])
    cb = fig2.colorbar(mappable=im3, cax=cax, ax=ax2)
    cb.set_label("$T_B$ (p-H$_2$CO $3_{0,3}-2_{0,2}$)")

    ax4 = pl.subplot(2,3,4)
    yi, xi = np.indices(mf.tline303.shape[1:])
    inds = [chi2b.argmin(axis=0), yi, xi]
    inds = [zz, slice(None), slice(None)]
    pl.imshow(mf.tline321[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              extent=mf.crange+mf.drange, vmin=vmin, vmax=vmax)
    pl.contour(mf.carr, mf.darr, chi2b.min(axis=0), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label_short)
    ax4.set_xticks(np.arange(11,16))
    pl.ylabel(density_label_short)

    ax5 = pl.subplot(2,3,5)
    zi, xi = np.indices([mf.tline303.shape[0], mf.tline303.shape[2],])
    inds = [zi, chi2b.argmin(axis=1), xi]
    inds = [slice(None), yy, slice(None)]
    pl.imshow(mf.tline321[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(mf.crange)/np.diff(mf.trange),
              extent=mf.crange+mf.trange, vmin=vmin, vmax=vmax)
    pl.contour(mf.carr, mf.tarr, chi2b.min(axis=1), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(column_label_short)
    ax5.set_xticks(np.arange(11,16))
    pl.ylabel(temperature_label)

    ax6 = pl.subplot(2,3,6)
    zi, yi = np.indices([mf.tline303.shape[0], mf.tline303.shape[2],])
    inds = [zi, yi, chi2b.argmin(axis=2)]
    inds = [slice(None), slice(None), xx]
    im6 = pl.imshow(mf.tline321[inds], cmap=pl.cm.bone_r, interpolation='spline36',
              norm=pl.matplotlib.colors.LogNorm(),
              aspect=np.diff(mf.drange)/np.diff(mf.trange),
              extent=mf.drange+mf.trange, vmin=vmin, vmax=vmax)
    pl.contour(mf.darr, mf.tarr,  chi2b.min(axis=2), levels=chi2b.min()+np.arange(nlevs))
    pl.xlabel(density_label_short)
    ax6.xaxis.set_ticks(np.arange(mf.darr.min(), mf.darr.max()))
    pl.ylabel(temperature_label)
    cax = fig2.add_axes([0.91,0.1,0.02,0.35])
    cb = fig2.colorbar(mappable=im6, cax=cax, ax=ax5)
    cb.set_label("$T_B$ (p-H$_2$CO $3_{2,1}-2_{2,0}$)")

    pl.suptitle(row['Source_Name'])
    pl.subplots_adjust(wspace=0.33, left=0.1, hspace=0.25)

    pl.savefig(paths.fpath('param_fits/{name}_{num}_h2coratio_minaxis.pdf'.format(name=row['Source_Name'],
                                                                                  num=num)), bbox_inches='tight')


    # shape is temperature, density, column
    # but we want the order shown here
    # so we're going to (very very very dangerously) swapaxes everything
    if hasattr(chi2r,'min'):
        chi2r = chi2r.swapaxes(1,2)
    if hasattr(chi2r2,'min'):
        chi2r2 = chi2r2.swapaxes(1,2)
    chi2b = chi2b.swapaxes(1,2)
    chi2X = chi2X.swapaxes(1,2)
    chi2_h2 = chi2_h2.swapaxes(1,2)
    chi2_ff = chi2_ff.swapaxes(1,2)
    chi2_dens = chi2_dens.swapaxes(1,2)
    tline303 = mf.tline303.swapaxes(1,2)

    for xax,yax,xlabel,ylabel,axis,ptype in zip((mf.darr,mf.darr,mf.carr),
                                                (mf.carr,mf.tarr,mf.tarr),
                                                (density_label, density_label,
                                                 column_label),
                                                (column_label,
                                                 temperature_label,
                                                 temperature_label),
                                                (0, 1, 2),
                                                ('dens_col','dens_tem','col_tem')):
        # Show the constraints provided by individual parameters
        fig3 = pl.figure(3)
        fig3.clf()
        # chi2b = chi2r + chi2_1 + chi2_2 + chi2X + chi2_h2
        ax1 = pl.subplot(2,2,1)
        if hasattr(chi2r, 'min'):
            pl.contourf(xax, yax, chi2r.min(axis=axis), levels=chi2r.min()+np.arange(nlevs), alpha=0.5)
        pl.contour(xax, yax, chi2b.min(axis=axis), levels=chi2b.min()+np.arange(nlevs))
        #if hasattr(chi2r2, 'min'):
        #    pl.contour(xax, yax, chi2r2.min(axis=axis),
        #               levels=chi2r2.min()+np.arange(nlevs),
        #               cmap=pl.cm.bone)
        pl.title(r"Ratio $R_1={0:0.3f}\pm{1:0.3f}$"
                 .format(ratio,eratio))
        ax4 = pl.subplot(2,2,2)
        pl.contourf(xax, yax, chi2X.min(axis=axis), levels=chi2X.min()+np.arange(nlevs), alpha=0.5)
        pl.contour(xax, yax, chi2b.min(axis=axis), levels=chi2b.min()+np.arange(nlevs))
        pl.title("log(p-H$_2$CO/H$_2$) $= {0:0.1f}\pm{1:0.1f}$".format(logabundance, elogabundance))
        ax3 = pl.subplot(2,2,3)
        #pl.contourf(xax, yax, chi2_h2.min(axis=axis), levels=chi2_h2.min()+np.arange(nlevs), alpha=0.5)
        pl.contourf(xax, yax, chi2_dens.min(axis=axis), levels=chi2_dens.min()+np.arange(nlevs), alpha=0.5)
        pl.contour(xax, yax, chi2b.min(axis=axis), levels=chi2b.min()+np.arange(nlevs))
        #pl.title("Total log$(N(\\mathrm{{H}}_2)) = {0:0.1f}\pm{1:0.1f}$".format(logh2column,
        #                                                                        elogh2column))
        pl.title("Density $n>10^{{{0:0.1f}}}$cm$^{{-3}}$".format(mindens))

        #ax5 = pl.subplot(2,2,4)
        #pl.contourf(xax, yax, (chi2_ff.min(axis=axis)),
        #            levels=chi2_ff.min()+np.arange(nlevs), alpha=0.5)
        #pl.contour(xax, yax, chi2b.min(axis=axis), levels=chi2b.min()+np.arange(nlevs))
        #pl.contour(xax, yax, (tline303 < 10*par1).max(axis=axis), levels=[0.5], colors='k')
        ##pl.contour(xax, yax, (tline303 < 100*par1).max(axis=axis), levels=[0.5], colors='k')
        ##pl.contour(xax, yax, (tline321 < 10*par2).max(axis=axis), levels=[0.5], colors='k', linestyles='--')
        ##pl.contour(xax, yax, (tline321 < 100*par2).max(axis=axis), levels=[0.5], colors='k', linestyles='--')
        #pl.xlabel(xlabel)
        #pl.ylabel(ylabel)
        #pl.title("Line Brightness + $ff\leq1$")

        ax5 = pl.subplot(2,2,4)
        ax5.contourf(xax, yax, chi2_h2.min(axis=axis), levels=chi2_h2.min()+np.arange(nlevs), alpha=0.5)
        ax5.contour(xax, yax, chi2b.min(axis=axis), levels=chi2b.min()+np.arange(nlevs))
        ax5.set_title("Total log$(N(\\mathrm{{H}}_2)) = {0:0.1f}\pm{1:0.1f}$"
                      .format(logh2column, elogh2column))

        pl.text(0.05, 0.5, ylabel, horizontalalignment='center',
                verticalalignment='center',
                rotation='vertical', transform=fig3.transFigure)
        pl.text(0.5, 0.02, xlabel, horizontalalignment='center', transform=fig3.transFigure)

        if xax is mf.carr:
            for ss in range(1,5):
                ax = pl.subplot(2,2,ss)
                ax.xaxis.set_ticks(np.arange(mf.carr.min(), mf.carr.max()))

        pl.subplots_adjust(wspace=0.2, hspace=0.4)
        outf = paths.fpath('param_fits/{name}_{ptype}_{num}_parameter_constraints.pdf'.format(name=row['Source_Name'],
                                                                                              ptype=ptype,
                                                                                              num=num))
        pl.savefig(outf, bbox_inches='tight')

    # IGNORE 321/322: it is generally not well constrained anyway
    mf.chi2 -= mf.chi2_r321322

    row_data = mf.get_parconstraints()
    for key,value in row_data.iteritems():
        row[key] = value


    width = row['width']*u.km/u.s
    row['reff_pc'] = reff.to(u.pc).value

    row['tkin_turb'] = heating.tkin_all(density=10**row['density_chi2']*u.cm**-3,
                                        sigma=width,
                                        lengthscale=reff,
                                        gradient=width/reff,
                                        tdust=row['higaldusttem']*u.K,
                                        crir=0./u.s)


    #if row_data['temperature_chi2'] == 10:
    #    import ipdb; ipdb.set_trace()

log.info("Completed source loop.")

fittable.write(os.path.join(analysispath,
                            'fitted_line_parameters_Chi2Constraints.ipac'),
               format='ascii.ipac')

log.info("Wrote table file.  Continuing to parameter plots.")


execfile(paths.pcpath('parameter_comparisons.py'))

pl.show()
