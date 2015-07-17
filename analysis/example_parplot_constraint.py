import pprint
import pylab as pl
import numpy as np
from constrain_parameters import paraH2COmodel
from astropy.utils.console import ProgressBar
import paths

mf = paraH2COmodel()

def make_model(density=4.5, column=13.5, temperature=75):
    dind = np.argmin(np.abs(mf.darr-density))
    cind = np.argmin(np.abs(mf.carr-column))
    tind = np.argmin(np.abs(mf.tarr-temperature))
    density = mf.darr[dind]
    column = mf.carr[cind]
    temperature = mf.tarr[tind]
    ta303 = mf.tline303[tind,dind,cind]
    ta321 = mf.tline321[tind,dind,cind]
    r321303 = ta321/ta303
    linewidth = 10
    logabundance = np.log10(1.2e-9)
    logh2column = column + np.log10(linewidth) + np.log10(1.) - logabundance 


    mf.set_constraints(ratio321303=r321303, eratio321303=0.01, logh2column=logh2column,
                       elogh2column=1, logabundance=logabundance,
                       elogabundance=1, mindens=4, linewidth=linewidth,
                       taline303=ta303, etaline303=0.07,
                       taline321=ta321, etaline321=0.07)
    constraints = mf.get_parconstraints()
    return constraints, mf, density, column, temperature

def make_model_plot(**kwargs):
    constraints,mf,density,column,temperature = make_model(**kwargs)


    pl.figure(1)
    mf.parplot1d('dens')
    pl.vlines(density, 0, pl.gca().get_ylim()[1], color='k', zorder=-10)
    pl.figure(2)
    mf.parplot1d('tem')
    pl.vlines(temperature, 0, pl.gca().get_ylim()[1], color='k', zorder=-10)
    pl.figure(3)
    mf.parplot1d('col')
    pl.vlines(column, 0, pl.gca().get_ylim()[1], color='k', zorder=-10)
    pl.figure(4)
    mf.parplot1d_all(levels=[0.68268949213708585])
    pl.subplot(3,1,1).vlines(column, 0, pl.gca().get_ylim()[1], color='k', zorder=-10)
    pl.subplot(3,1,2).vlines(density, 0, pl.gca().get_ylim()[1], color='k', zorder=-10)
    pl.subplot(3,1,3).vlines(temperature, 0, pl.gca().get_ylim()[1], color='k', zorder=-10)
    pl.suptitle('t={0:0.1f}_c={1:0.1f}_n={2:0.1f}.png'.format(temperature,
                                                              column,
                                                              density,))

    pl.savefig(paths.fpath('param_fits/oned_fit_parameters_example_t={0:0.1f}_c={1:0.1f}_n={2:0.1f}.png'.format(temperature,
                                                                                                                column,
                                                                                                                density,)),
               bbox_inches='tight')
    pl.figure(5)
    mf.denstemplot()
    pl.plot(density, temperature, 'o', markerfacecolor='none', markeredgecolor='r', markeredgewidth=2)
    pl.figure(7)
    mf.denscolplot()
    pl.plot(density, column, 'o', markerfacecolor='none', markeredgecolor='r', markeredgewidth=2)
    pl.figure(8)
    mf.coltemplot()
    pl.plot(column, temperature, 'o', markerfacecolor='none', markeredgecolor='r', markeredgewidth=2)

    pl.draw()
    pl.show()
    return constraints,mf


def make_many_models():

    constraints = {}
    for column in (12.5, 13, 13.5, 14):
        for density in (4,4.5,5,5.5):
            for temperature in (25,50,75,100):
                for ii in range(1,9): pl.figure(ii).clf()
                fits,mf = make_model_plot(density=density, column=column, temperature=temperature)
                constraints[(column, density, temperature)] = fits

    return constraints

def diffplot_t_of_n(densities=np.linspace(3,7), column=13.5, temperatures=[25,50,75,100,125,150]):
    grid_exp = np.empty([len(densities), len(temperatures)])
    grid_ML = np.empty([len(densities), len(temperatures)])
    ngrid_exp = np.empty([len(densities), len(temperatures)])
    ngrid_ML = np.empty([len(densities), len(temperatures)])
    for idens, density in enumerate(ProgressBar(densities)):
        for item, temperature in enumerate(temperatures):
            constraints,mf,density,column,temperature = make_model(density=density, temperature=temperature, column=column)
            grid_exp[idens,item] = constraints['expected_temperature']
            grid_ML[idens,item] = constraints['temperature_chi2']
            ngrid_exp[idens,item] = constraints['expected_density']
            ngrid_ML[idens,item] = constraints['density_chi2']

    pl.figure(1).clf()

    for ii,(tem,color) in enumerate(zip(temperatures,('r','g','b','c','m','orange'))):
        pl.plot(densities, grid_exp[:,ii], color=color)
        pl.plot(densities, grid_ML[:,ii], '--', color=color)
        pl.hlines(tem, densities.min(), densities.max(), label='T={0}K'.format(tem), color=color)
    pl.plot([], 'k', label='Expectation Value')
    pl.plot([], 'k--', label='Maximum Likelihood')
    pl.xlabel("log $n$(H$_2$) [cm$^{-3}$]")
    pl.ylabel("Temperature (K)")
    pl.legend(loc='best', fontsize=14)

    pl.figure(2).clf()

    for ii,(tem,color) in enumerate(zip(temperatures,('r','g','b','c','m','orange'))):
        pl.plot(densities, (grid_exp[:,ii]-tem)/tem, color=color, label='T={0}K'.format(tem))
        pl.plot(densities, (grid_ML[:,ii]-tem)/tem, '--', color=color)
    pl.plot([], 'k', label='Expectation Value')
    pl.plot([], 'k--', label='Maximum Likelihood')
    pl.xlabel("log n(H$_2$CO) [cm$^{-3}$]")
    pl.ylabel("Fractional Difference\n(recovered-input)/input")
    pl.legend(loc='best', fontsize=14)
    pl.ylim(-0.5,0.5)
    pl.grid()


    pl.figure(3).clf()

    for ii,(tem,color) in enumerate(zip(temperatures,('r','g','b','c','m','orange'))):
        pl.plot(densities, ngrid_exp[:,ii], color=color)
        pl.plot(densities, ngrid_ML[:,ii], '--', color=color)
        pl.plot(densities, densities, label='T={0}K'.format(tem), color=color)
    pl.plot([], 'k', label='Expectation Value')
    pl.plot([], 'k--', label='Maximum Likelihood')
    pl.xlabel("Input log $n$(H$_2$) [cm$^{-3}$]")
    pl.ylabel("Recovered log $n$(H$_2$) [cm$^{-3}$]")
    pl.legend(loc='best', fontsize=14)

    pl.figure(4).clf()

    for ii,(tem,color) in enumerate(zip(temperatures,('r','g','b','c','m','orange'))):
        pl.plot(densities, (10**ngrid_exp[:,ii]-10**densities)/10**densities, color=color, label='T={0}K'.format(tem))
        pl.plot(densities, (10**ngrid_ML[:,ii] -10**densities)/10**densities, '--', color=color)
    pl.plot([], 'k', label='Expectation Value')
    pl.plot([], 'k--', label='Maximum Likelihood')
    pl.xlabel("Input log $n$(H$_2$) [cm$^{-3}$]")
    pl.ylabel("Fractional Difference\n(recovered-input)/input")
    pl.legend(loc='best', fontsize=14)
    pl.ylim(-0.5,0.5)
    pl.grid()

    return densities, grid_exp, grid_ML, ngrid_exp, ngrid_ML

def diffplot_t_of_c(density=4.5, columns=np.linspace(12,15), temperatures=[25,50,75,100,125,150]):
    grid_exp = np.empty([len(columns), len(temperatures)])
    grid_ML = np.empty([len(columns), len(temperatures)])
    for icol, column in enumerate(ProgressBar(columns)):
        for item, temperature in enumerate(temperatures):
            constraints,mf,density,column,temperature = make_model(density=density, temperature=temperature, column=column)
            grid_exp[icol,item] = constraints['expected_temperature']
            grid_ML[icol,item] = constraints['temperature_chi2']

    pl.figure(1).clf()

    for ii,(tem,color) in enumerate(zip(temperatures,('r','g','b','c','m','orange'))):
        pl.plot(columns, grid_exp[:,ii], color=color)
        pl.plot(columns, grid_ML[:,ii], '--', color=color)
        pl.hlines(tem, columns.min(), columns.max(), label='T={0}K'.format(tem), color=color)
    pl.plot([], 'k', label='Expectation Value')
    pl.plot([], 'k--', label='Maximum Likelihood')
    pl.xlabel("log N(H$_2$CO) [cm$^{-2}$]")
    pl.ylabel("Temperature (K)")
    pl.legend(loc='best', fontsize=14)

    pl.figure(2).clf()

    for ii,(tem,color) in enumerate(zip(temperatures,('r','g','b','c','m','orange'))):
        pl.plot(columns, (grid_exp[:,ii]-tem)/tem, color=color, label='T={0}K'.format(tem))
        pl.plot(columns, (grid_ML[:,ii]-tem)/tem, '--', color=color)
    pl.plot([], 'k', label='Expectation Value')
    pl.plot([], 'k--', label='Maximum Likelihood')
    pl.xlabel("log N(H$_2$CO) [cm$^{-2}$]")
    pl.ylabel("Fractional Difference\n(recovered-input)/input")
    pl.legend(loc='best', fontsize=14)
    pl.ylim(-0.5,0.5)
    pl.grid()

    return columns, grid_exp, grid_ML


def main():
    make_model()
    
    densities, dgrid_exp, dgrid_ML, dngrid_exp, dngrid_ML = diffplot_t_of_n()
    pl.figure(1).savefig(paths.fpath('param_fits/modeled_temperature_vs_density_real_vs_recovered_C=13.5.png'), bbox_inches='tight')
    pl.figure(2).savefig(paths.fpath('param_fits/modeled_temperature_vs_density_real_vs_recovered_fractions_C=13.5.png'), bbox_inches='tight')
    pl.figure(3).savefig(paths.fpath('param_fits/modeled_density_vs_density_real_vs_recovered_C=13.5.png'), bbox_inches='tight')
    pl.figure(4).savefig(paths.fpath('param_fits/modeled_density_vs_density_real_vs_recovered_fractions_C=13.5.png'), bbox_inches='tight')

    columns, cgrid_exp, cgrid_ML = diffplot_t_of_c()
    pl.figure(1).savefig(paths.fpath('param_fits/modeled_temperature_vs_column_real_vs_recovered_n=4.5.png'), bbox_inches='tight')
    pl.figure(2).savefig(paths.fpath('param_fits/modeled_temperature_vs_column_real_vs_recovered_fractions_n=4.5.png'), bbox_inches='tight')

    return locals()
