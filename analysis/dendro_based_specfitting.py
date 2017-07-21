"""
Use dendrograms as seeds...

This works, but it's not 100% reliable, which means it's going to be tough to
analyze.  Also, it seems that there are some spectra that don't really want to
fit.  You give 'em guesses, and they just sit there.  I've seen some really
strange baseline solutions that are pretty much inexplicable (n'th order = 1st
order always?).

There may also be a problem with the parallelization; it is not clear that the
fits are independent.

"""
import os
from astropy.table import Table
from pyspeckit_fitting import simple_fitter, simple_fitter2
from pyspeckit.parallel_map import parallel_map
from full_cubes import (pcube_merge_high, cube_merge_high, pcube_merge_high_sm,
                        cube_merge_high_sm)
from paths import mpath,hpath,fpath
from dendrograms import dend,dendsm,catalog,catalog_sm
import numpy as np
from astropy.io import fits
from astropy.utils.console import ProgressBar
from astropy import log
#from agpy.mad import MAD
from astropy.stats import mad_std as MAD
import multiprocessing
import pyspeckit

def unique_indices(xyinds):
    raise NotImplementedError("This approach fails because it sorts *BOTH* axes independently")
    assert xyinds.shape[1] == 2
    xyinds.sort(axis=0)
    d = np.diff(xyinds, axis=0)
    ok = np.sum(d, axis=1) > 0
    inds = np.concatenate([[xyinds[0]], xyinds[1:,:][ok,:]])
    return inds

def get_all_indices(dendrogram):
    inds = []
    for structure in dendrogram.trunk:
        if structure.parent is None: # should be true for all
            inds.append(structure.indices())
    allzinds,allyinds,allxinds = [np.concatenate(ii)
                                  for ii in zip(*inds)]
    unique_inds = set(zip(allyinds, allxinds))
    return np.array(list(unique_inds))

def match_position(position, structure, **kwargs):
    if hasattr(position,'ndim') and position.ndim == 2:
        return overlaps_positions(position, structure, **kwargs)
    y,x = position
    #return (x,y) in zip(*structure.indices()[1:])
    #return (x in structure.indices()[2] and y in structure.indices()[1])
    inds = structure.indices()
    return np.count_nonzero((x == inds[2]) & (y == inds[1])) > 0

def overlaps_positions(positions, structure, pxmin=None, pxmax=None,
                       pymin=None, pymax=None):
    if pymin is None or pxmin is None:
        pymin,pxmin = np.array(positions).min(axis=1)
    if pymax is None or pxmax is None:
        pymax,pxmax = np.array(positions).max(axis=1)

    inds = structure.indices()
    xmin,xmax = inds[2].min(),inds[2].max()
    ymin,ymax = inds[1].min(),inds[1].max()

    if pxmax < xmin:
        log.debug("pxmax < xmin: {0},{1}".format(pxmax,xmin))
        return False
    if pxmin > xmax:
        log.debug("pxmin > xmax: {0},{1}".format(pxmin,xmax))
        return False
    if pymax < ymin:
        log.debug("pymax < ymin: {0},{1}".format(pymax,ymin))
        return False
    if pymin > ymax:
        log.debug("pymin > ymax: {0},{1}".format(pymin,ymax))
        return False

    for y,x in positions.T:
        if np.count_nonzero((x == inds[2]) & (y == inds[1])) > 0:
            return True

    return False



def smallest_match(position, structure, **kwargs):
    """
    Given a position and a structure in which at least one structure in the
    tree contains a match, return the smallest structure that matches that
    spatial position
    """
    print(structure)
    if match_position(position, structure, **kwargs):
        if structure.is_leaf:
            return structure
        else:
            c1 = smallest_match(position, structure.children[0])
            if c1 is not structure:
                return c1
            return smallest_match(position, structure.children[1])
    else:
        return structure.parent

def all_matches(position, structure, **kwargs):
    """
    Given a position and a structure in which at least one structure in the
    tree contains a match, return all structures that overlap with that
    spatial position
    """
    if match_position(position, structure, **kwargs):
        if structure.children:
            return all_matches(position, structure.children[0]) +\
                   all_matches(position, structure.children[1]) +\
                   [structure]
        else:
            return [structure]
    else:
        return []

def excise_parents(structure_list):
    to_remove = []
    for obj in structure_list:
        while obj.parent in structure_list:
            to_remove.append(obj.parent)
            obj = obj.parent
    return [obj for obj in structure_list if obj not in to_remove]

def excise_parents_of_object(structure_list, obj):
    to_remove = []
    while obj.parent in structure_list:
        to_remove.append(obj.parent)
        obj = obj.parent
    return [obj for obj in structure_list if obj not in to_remove]

def excise_children(structure_list):
    to_remove = []
    for obj in structure_list:
        if obj.parent in structure_list:
            to_remove.append(obj)
    return [obj for obj in structure_list if obj not in to_remove]

def find_overlapping_trunks(position, dendrogram, **kwargs):
    return [d for d in dendrogram.trunk
            if d.parent is None and match_position(position, d, **kwargs)]

def find_smallest_overlapping_branches(position, dendrogram):
    if hasattr(position,'ndim') and position.ndim == 2:
        pymin,pxmin = np.array(position).min(axis=1)
        pymax,pxmax = np.array(position).max(axis=1)
    else:
        pymin,pymax = position[0],position[0]
        pxmin,pxmax = position[0],position[0]
    matched_roots = find_overlapping_trunks(position, dendrogram, pxmin=pxmin,
                                            pxmax=pxmax, pymin=pymin,
                                            pymax=pymax)
    matched_branches = [smallest_match(position, d, pxmin=pxmin, pxmax=pxmax,
                                       pymin=pymin, pymax=pymax)
                        for d in matched_roots]
    return matched_branches

def find_all_overlapping_branches(position, dendrogram):
    matched_roots = find_overlapping_trunks(position, dendrogram)
    matched_branches = [a
                        for d in matched_roots
                        for a in all_matches(position, d) 
                       ]
    return matched_branches

def get_guesses(index, catalog, max_comp=3):
    """
    put the guesses into flat format and limit their numbers
    """
    amp,vc,vsig,r = catalog['Smean303', 'v_cen', 'v_rms',
                            'r303321'][index].columns.values()
    keep_velo = (vsig < 10) & (vsig > 1)
    if np.count_nonzero(keep_velo) > max_comp:
        # Also exclude small objects when there are many overlaps
        big_objs = catalog[index]['npix'] > 100
        # but only if excluding them isn't overly restrictive
        if np.count_nonzero(big_objs & keep_velo) >= max_comp:
            keep_velo &= big_objs
    keep_inds = np.argsort(amp[keep_velo])[-max_comp:]
    log.debug('Kept {0}, amps {1}'.format(np.array(index)[keep_velo][keep_inds],
                                          np.array(amp)[keep_velo][keep_inds]))
    amp,vc,vsig,r = [a[keep_velo][keep_inds] for a in (amp,vc,vsig,r)]
    glon = [xc - 360 if xc > 180 else xc
            for xc in catalog['x_cen']]
    # TODO: make sure the velocity units remain consistent
    # The logic here is to filter out things at vlsr<-50 km/s that are not in Sgr C;
    # these are mostly other lines in Sgr B2
    result = [pars
              for pars,glon in zip(zip(*[x.data for x in [amp,vc/1e3,vsig,r,amp]]),
                                      glon)
             ]

    return result

def clean_catalog(catalog):
    """
    Remove cataloged clumps that are not H2CO
    """
    bad = ((catalog['v_cen'] < -50) &
           (catalog['x_cen']>0) &
           (catalog['x_cen']<180))
    raise # this is a bad idea (makes catalog/dend disagree)
    return catalog[~bad]

pcube_merge_high.Registry.add_fitter('h2co_simple', simple_fitter, 5,
                                     multisingle='multi')
pcube_merge_high.Registry.add_fitter('h2co_simple2', simple_fitter2, 6,
                                     multisingle='multi')
pcube_merge_high_sm.Registry.add_fitter('h2co_simple', simple_fitter, 5,
                                     multisingle='multi')
pcube_merge_high_sm.Registry.add_fitter('h2co_simple2', simple_fitter2, 6,
                                     multisingle='multi')

def fit_position(position, dendrogram=dend, pcube=pcube_merge_high,
                 catalog=catalog, **kwargs):
    branches = find_smallest_overlapping_branches(position, dendrogram)
    branch_ids = [l.idx for l in branches]
    if len(branch_ids) == 0:
        log.error("Invalid position given: {0}.  No branches found.".format(position))
        raise
        return
    guess = get_guesses(branch_ids, catalog)

    sp = pcube.get_spectrum(position[1], position[0])
    sp.specname = '{1},{0}'.format(*position)

    return fit_spectrum(sp, guess, **kwargs)

def fit_object(obj, dendrogram=dend, pcube=pcube_merge_high, catalog=catalog,
               **kwargs):

    xyinds = np.array(obj.indices()[1:])
    mean_position = xyinds.mean(axis=1).astype('int')
    u_positions = np.array(list(set(zip(*xyinds))))
    #if mean_position not in u_positions:
    #    raise NotImplementedError("Find the nearest pixel that IS in the inds...")
    #assert match_position(mean_position, obj)

    log.debug("Locating overlapping branches....")
    all_branches = find_all_overlapping_branches(u_positions.T, dendrogram)
    branches = excise_children(excise_parents_of_object(all_branches, obj))
    branch_ids = [l.idx for l in branches if not l.bad]
    if len(branch_ids) == 0:
        log.error("Invalid position given: {0}.  No branches"
                  " found.".format(mean_position))
        raise
        return
    log.info("Found overlapping branches {0}".format(branch_ids))
    guess = get_guesses(branch_ids, catalog)
    log.info("Guess: {0}".format(guess))
    if len(guess) > 20:
        raise ValueError("Too many guesses.  Too many components.")

    # This needs replacing...
    # we want to extract the full spectrum...
    dend_obj_mask = obj.get_mask()
    all_pix = dend_obj_mask.max(axis=0)
    ia,ib = np.where(all_pix)
    view = slice(None), slice(ia.min(), ia.max()), slice(ib.min(), ib.max())
    weight = dend_obj_mask[view].sum(axis=0)
    assert weight.sum() > 0
    sp_data = np.nansum(pcube_merge_high.cube[view] * weight, axis=(1,2)) / weight.sum()

    sp = pyspeckit.Spectrum(data=sp_data, xarr=pcube_merge_high.xarr.as_unit('GHz'),
                            error=np.zeros_like(sp_data)+MAD(sp_data),
                            header=fits.Header())
    sp.specname = "ID{0}_{1},{2}".format(obj.idx, mean_position[1], mean_position[0])
    sp.Registry.add_fitter('h2co_simple', simple_fitter, 5,
                           multisingle='multi')
    sp.Registry.add_fitter('h2co_simple2', simple_fitter2, 6,
                           multisingle='multi')

    log.debug("Spectrum loaded, now fitting.")

    return fit_spectrum(sp, guess, **kwargs)

def fit_spectrum(sp, guess, plot=True, order=1, second_ratio=False,
                 fig=pl.figure(1),
                 verbose=False, **kwargs):
    """
    Quasi-automatic fit to a spectrum...
    """

    if second_ratio:
        guesses = [p
                   for g in guess if g[0]>0 and g[3]>0
                   for p in (g[:4]+(1,g[4]))]
        limits = [a
                  for g in guess if g[0]>0 and g[3]>0
                  for a in ((g[0]*0.1, g[0]/0.1),
                            (g[1]-15, g[1]+15),
                            (max(g[2]-5,1), g[2]+10),
                            (g[3]*0.1, g[3]/0.1),
                            (0.2, 1.4), # physical limits (0.3,1.1)
                            (g[4]*0.05, g[4]/0.05))
                 ]
        fittype = 'h2co_simple2'
        assert len(guesses) % 6 == 0
        assert len(limits) % 6 == 0
        if len(guesses) > 0:
            assert guesses[0] in [g[0] for g in guess]
        n = 6
    else:
        guesses = [p for g in guess if g[0]>0 and g[3]>0 for p in g]
        limits = [a
                  for g in guess if g[0]>0 and g[3]>0
                  for a in ((g[0]*0.1, g[0]/0.1),
                            (g[1]-15, g[1]+15),
                            (max(g[2]-5,1), g[2]+10),
                            (g[3]*0.1, g[3]/0.1),
                            (g[4]*0.05, g[4]/0.05))
                 ]
        fittype = 'h2co_simple'
        assert len(guesses) % 5 == 0
        assert len(limits) % 5 == 0
        n = 5

    # Widths could be initialized out of range
    for ii in range(2,len(guesses),n):
        if guesses[ii] < limits[ii][0]:
            guesses[ii] = 1.1
        if guesses[ii] > limits[ii][1]:
            guesses[ii] = limits[ii][1]*0.9

    if len(guesses) == 0:
        #log.info("Position {0},{1} has no guesses.".format(position[0],
        #                                                   position[1]))
        log.debug("Guess was: {0}".format(guess))
        return

    # Some positions have bad noise values; these are not worth wasting time on
    if all(np.isnan(sp.error)):
        return

    sp.specfit(fittype=fittype, multifit=True,
               guesses=guesses,
               limited=[(True,True)] * len(guesses),
               limits=limits,
               verbose=verbose,
               renormalize=False,
               **kwargs
              )

    assert len(guesses) % (5+second_ratio) == 0
    assert len(sp.specfit.parinfo) % (5+second_ratio) == 0

    if plot:
        sp.plotter(figure=fig)
        sp.specfit.plot_fit()
        sp.baseline(excludefit=True, subtract=True, highlight_fitregion=True, order=order)
    else:
        sp.baseline(excludefit=True, subtract=True, order=order, verbose=verbose)
        
    sp.specfit(fittype=fittype, multifit=True,
               guesses=guesses,
               limited=[(True,True)] * len(guesses),
               limits=limits,
               verbose=verbose,
               renormalize=False,
               **kwargs
              )

    assert len(guesses) % (5+second_ratio) == 0
    assert len(sp.specfit.parinfo) % (5+second_ratio) == 0

    if plot:
        sp.plotter(figure=fig)
        sp.specfit.plot_fit(show_components=True)
        #sp.baseline.plot_baseline()

    return sp

def fit_all_objects(dend=dend, **kwargs):
    fits = {}

    for obj in dend:
        if obj.bad:
            log.info("Skipping object {0} because it is HC3N or other.".format(obj.idx))
            continue
        log.info("Operating on object {0}".format(obj.idx))
        try:
            sp = fit_object(obj, plot=False, **kwargs)
        except ValueError as ex:
            log.info("Exception for {0}: {1}".format(obj.idx, ex.message))
            continue
        if sp is None:
            log.info("No result for {0}".format(obj.idx))
            continue
        sp.plotter(figure=pl.figure(1))
        sp.specfit.plot_fit(show_components=True)
        sp.plotter.savefig(fpath('dendro/dendro_fit_obj{0:04d}.png'.format(obj.idx)))
        fits[obj.idx] = sp.specfit.parinfo

    return fits

def fit_all_positions(dendrogram=dend, pcube=pcube_merge_high, catalog=catalog,
                      order=1, second_ratio=False, ncores=1, positions=None,
                      outfilename=None):
    if positions is None:
        # Reverse order: start from the smallest trees
        # Positions are y,x
        positions = get_all_indices(dendrogram)[::-1]

    log.info("Fitting {0} positions.".format(len(positions)))

    if outfilename is not None:
        fitted_positions,parvalues,parerrors = read_pars(outfilename)
        outfile = True
    else:
        fitted_positions = []
        outfile = None

    lock = multiprocessing.Lock()

    def get_fitvals(p, plot=False, order=order, second_ratio=second_ratio,
                    outfile=outfile, lock=lock):
        if tuple(p) in fitted_positions:
            return

        log.debug("Fitting position {0}".format(p))
        result = fit_position(p, dendrogram=dendrogram, catalog=catalog,
                              pcube=pcube,
                              plot=False, order=order,
                              second_ratio=second_ratio)

        fitted_positions.append(tuple(p))
        if result is None:
            parvalues.append(None)
            parerrors.append(None)
            if outfile is not None:
                with lock:
                    with open(outfilename, 'a') as outfile:
                        outfile.write("{0}, {1}, {2}, {3}\n".format(p[0], p[1], None, None))
                        outfile.flush()
            return
        else:
            parvalues.append(result.specfit.parinfo.values)
            parerrors.append(result.specfit.parinfo.errors)
            if outfile is not None:
                with lock:
                    with open(outfilename, 'a') as outfile:
                        outfile.write("{0}, {1}, {2}, {3}\n".format(p[0], p[1],
                                                                    result.specfit.parinfo.values,
                                                                    result.specfit.parinfo.errors))
                        outfile.flush()
            return result.specfit.parinfo.values, result.specfit.parinfo.errors

    if ncores == 1:
        results = [get_fitvals(p, plot=False, order=order,
                               second_ratio=second_ratio)
                   for p in ProgressBar(positions)]
    else:
        results = parallel_map(get_fitvals, positions, numcores=ncores)

    bad_positions = [p for p,r in zip(positions,results) if r is None]
    positions2 = [p for p,r in zip(positions,results) if r is not None]
    results2 = [r for r in results if r is not None]

    return positions2,results2,bad_positions

def pars_to_maps(positions, parvalues, shape=pcube_merge_high.cube.shape[1:],
                 celwcs=cube_merge_high.wcs.celestial, suffix=""):

    maxlen = max(len(r) for r in parvalues if r is not None)

    if maxlen % 6 == 0:
        names = ['AMPLITUDE', 'VELOCITY', 'WIDTH', 'RATIO321303X',
                 'RATIO322321X', 'AMPCH3OH']
        n = 6
    elif maxlen % 5 == 0:
        names = ['AMPLITUDE', 'VELOCITY', 'WIDTH', 'RATIO321303X',
                 'AMPCH3OH']
        n = 5
    else:
        raise

    maps = []
    for jj in range(maxlen / n):
        for ii in xrange(n):
            arr = np.zeros(shape)
            for p,r in zip(positions,parvalues):
                ind = ii + jj*n
                if r is None:
                    continue
                elif len(r) > ind:
                    arr[p[1], p[0]] = r[ind]
            maps.append(arr)
            hdu = fits.PrimaryHDU(data=arr,
                                  header=celwcs.to_header())
            hdu.writeto(hpath("pyspeckit_{0}{1}{2}.fits".format(names[ii],
                                                                jj,
                                                                suffix)),
                        clobber=True)

    return maps

def read_pars(filename):
    fitted_positions,parvalues,parerrors = [],[],[]
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            for line in f.readlines():
                x,y,a,b = eval(line.strip())
                fitted_positions.append((x,y))
                parvalues.append(a)
                parerrors.append(b)
    return fitted_positions, parvalues, parerrors

def do_fitting(ncores=4):

    # smooth both
    results3 = fit_all_positions(dendrogram=dendsm, catalog=catalog_sm,
                                            pcube=pcube_merge_high_sm,
                                            second_ratio=True,
                                            outfilename=hpath('pyspeckit_fits_smsm.txt'),
                                            ncores=ncores)
    (positions_sm2, results_sm2,
     bad_positions_sm2) = read_pars(hpath('pyspeckit_fits_smsm.txt'))
    pars_to_maps(positions_sm2, results_sm2, suffix='_sm2')

    # sharp both
    results2 = fit_all_positions(dendrogram=dend, catalog=catalog,
                                        second_ratio=True,
                                        outfilename=hpath('pyspeckit_fits.txt'),
                                        ncores=ncores)
    (positions, results,
     bad_positions) = read_pars(hpath('pyspeckit_fits.txt'))
    pars_to_maps(positions, results, suffix='')

    # Smooth dendrograms, sharp image
    results = fit_all_positions(dendrogram=dendsm, catalog=catalog_sm,
                                second_ratio=True,
                                outfilename=hpath('pyspeckit_fits_densm.txt'),
                                ncores=ncores)
    (positions_sm1, results_sm1,
     bad_positions_sm1) = read_pars(hpath('pyspeckit_fits_densm.txt'))
    pars_to_maps(positions_sm1, results_sm1, suffix='_sm1')

