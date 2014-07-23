from shfi_otf_pipeline.make_apex_cubes import efuncs
import numpy as np
import pylab as pl

x = np.random.randn(200,100)

efuncarr_sp_all, covmat_sp_all, evals_sp_all, evects_sp_all = efuncs(x, return_others=True, neig=x.shape[1])

efuncarr_np_all, covmat_np_all, evals_np_all, evects_np_all = efuncs(x, return_others=True, neig=None)

efuncarr_sp_3, covmat_sp_3, evals_sp_3, evects_sp_3 = efuncs(x, return_others=True, neig=3)

efuncarr_sp_all[:, 3:] = 0
efuncarr_np_all[:, 3:] = 0

to_subtract_np_all = np.inner(efuncarr_np_all,evects_np_all).T
to_subtract_sp_all = np.inner(efuncarr_sp_all,evects_sp_all).T
to_subtract_sp_3 = np.inner(efuncarr_sp_3,evects_sp_3).T

pl.figure(1)
pl.clf()
pl.subplot(3,1,1); pl.imshow(to_subtract_np_all); pl.colorbar()
pl.subplot(3,1,2); pl.imshow(to_subtract_sp_all); pl.colorbar()
pl.subplot(3,1,3); pl.imshow(to_subtract_sp_3); pl.colorbar()

pl.figure(2)
pl.clf()
pl.subplot(3,1,1); pl.imshow(to_subtract_np_all-to_subtract_sp_3); pl.colorbar()
pl.subplot(3,1,2); pl.imshow(to_subtract_sp_all-to_subtract_sp_3); pl.colorbar()
pl.subplot(3,1,3); pl.imshow(to_subtract_np_all-to_subtract_sp_all); pl.colorbar()

pl.figure(3)
pl.clf()
pl.plot(evals_sp_all, label='spall')
pl.plot(evals_np_all, label='npall')
pl.plot(evals_sp_3, label='sp3')

pl.figure(4)
pl.clf()
for ii in range(3):
    pl.subplot(3,1,ii+1)
    pl.plot(evects_sp_all[:,ii])
    pl.plot(evects_np_all[:,ii])
    pl.plot(evects_sp_3[:,ii])


for covmat,w,v in [(covmat_np_all, evals_np_all, evects_np_all),
                   (covmat_sp_all, evals_sp_all, evects_sp_all),
                   (covmat_sp_3,   evals_sp_3, evects_sp_3)]:
    print (np.dot(covmat, v[:, 0]) - w[0] * v[:, 0]).sum()
