import pylab as pl
from astropy.io import fits
import paths
import numpy as np

ftemplate = '{line}_cube_smooth_bl_{ttype}.fits'

d322 = {'ttype':'wtdmeantemperature',
        'line':'H2CO_322221_to_303202'}
t322 = fits.getdata(paths.dpath(ftemplate.format(**d322)))

d321 = {'ttype':'wtdmeantemperature',
        'line':'H2CO_321220_to_303202'}
t321 = fits.getdata(paths.dpath(ftemplate.format(**d321)))

fig = pl.figure(1)
fig.clf()

histo321 = np.sort(t321, axis=0)
ax = fig.add_subplot(2,1,1)
im = ax.imshow(histo321, cmap='RdYlBu_r', vmin=20, vmax=120)
fig.colorbar(im)

histo322 = np.sort(t322, axis=0)
ax = fig.add_subplot(2,1,2)
im = ax.imshow(histo322, cmap='RdYlBu_r', vmin=20, vmax=120)
fig.colorbar(im)

fig2 = pl.figure(2)
#for threshold, sign, color in zip([20,30,40,50,50,80,100,120],
#                                 [-1,-1,-1,-1, 1, 1,  1,  1],
#                                 pl.cm.RdYlBu_r(np.linspace(0,1,8))):
#    pl.plot((t322*sign > threshold*sign).sum(axis=0), color=color)
for threshold, color in zip([20,30,40,50,50,80,100,120],
                            pl.cm.RdYlBu_r(np.linspace(0,1,8))):
    pl.plot((t322 > threshold).sum(axis=0), color=color)

pl.draw()
pl.show()
