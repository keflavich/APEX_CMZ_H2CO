from shfi_otf_pipeline import make_apex_cubes
import time # for performance analysis
from astropy import units as u
import FITS_tools
import image_tools
import numpy as np
# for bigger figures
import pylab as pl
import spectral_cube
from astropy.io import fits
import aplpy
import os
pl.rcParams['figure.figsize'] = (12.0, 10.0)

remake_cubes = False
if remake_cubes:

    make_apex_cubes.build_cube_2013(datasets=['M-091.F-0019-2013-2013-06-11',],
                                    lowhigh='high',
                                    pca_clean=True,
                                    pcakwargs={'smoothing_scale':25.,
                                               'ncomponents':3},
                                    extra_suffix='_3PCA')

    make_apex_cubes.build_cube_2013(datasets=['M-091.F-0019-2013-2013-06-11',],
                                    lowhigh='high',
                                    pca_clean=False,
                                    extra_suffix='_noPCA')

    make_apex_cubes.build_cube_2013(datasets=['M-091.F-0019-2013-2013-06-11',],
                                    lowhigh='high',
                                    pca_clean=True,
                                    pcakwargs={'smoothing_scale':25.,                                          
                                               'ncomponents':1},                                
                                    extra_suffix='_1PCA')

v1,v2 = 6*u.km/u.s, 114*u.km/u.s

dpath = make_apex_cubes.june2013path
cube = spectral_cube.SpectralCube.read(os.path.join(dpath,'APEX_H2CO_2013_high_noPCA_sub.fits'))
h2cocube = cube.with_spectral_unit(u.km/u.s, rest_value=218.22219*u.GHz, velocity_convention='radio')
h2cocubecut = h2cocube.spectral_slab(v1,v2)
inth2co = h2cocubecut.moment0()
peakh2co = fits.PrimaryHDU(h2cocubecut.max(axis=0).value, header=inth2co.hdu.header)


# In[8]:

F = aplpy.FITSFigure(inth2co.hdu)
F._ax1.set_title("Integrated h2co")
F.show_grayscale()
F.add_colorbar()

F = aplpy.FITSFigure(peakh2co)
F._ax1.set_title("Peak h2co")
F.show_grayscale()
F.add_colorbar()


# In[9]:

pcacube = spectral_cube.SpectralCube.read(os.path.join(dpath,'APEX_H2CO_2013_high_3PCA_sub.fits'))
h2copcacube = pcacube.with_spectral_unit(u.km/u.s, rest_value=218.22219*u.GHz, velocity_convention='radio')
h2copcacubecut = h2copcacube.spectral_slab(v1, v2)
pcainth2co = h2copcacubecut.moment0()
pcapeakh2co = fits.PrimaryHDU(h2copcacubecut.max(axis=0).value, header=pcainth2co.hdu.header)


# In[10]:

F = aplpy.FITSFigure(pcainth2co.hdu)
F._ax1.set_title("Integrated h2co")
F.show_grayscale()
F.add_colorbar()

F = aplpy.FITSFigure(pcapeakh2co)
F._ax1.set_title("Peak h2co")
F.show_grayscale()
F.add_colorbar()


# In[11]:

pca1cube = spectral_cube.SpectralCube.read(os.path.join(dpath,'APEX_H2CO_2013_high_1PCA_sub.fits'))
h2copca1cube = pca1cube.with_spectral_unit(u.km/u.s, rest_value=218.22219*u.GHz, velocity_convention='radio')
h2copca1cubecut = h2copca1cube.spectral_slab(v1, v2)
pca1inth2co = h2copca1cubecut.moment0()
pca1peakh2co = fits.PrimaryHDU(h2copca1cubecut.max(axis=0).value, header=pca1inth2co.hdu.header)


# In[12]:

F = aplpy.FITSFigure(pca1inth2co.hdu)
F._ax1.set_title("Integrated h2co")
F.show_grayscale()
F.add_colorbar()

F = aplpy.FITSFigure(pca1peakh2co)
F._ax1.set_title("Peak h2co")
F.show_grayscale()
F.add_colorbar()


# In[13]:

diffh2co = fits.PrimaryHDU(inth2co.hdu.data - pcainth2co.hdu.data, header=inth2co.hdu.header)
diffh2co2 = fits.PrimaryHDU(inth2co.hdu.data - pca1inth2co.hdu.data, header=inth2co.hdu.header)


# In[14]:

fig = pl.figure(figsize=(16,8))
subplot_diff2 = 0.4
subplot_bounds2 = [0.05,0.5]
subplot_diff3 = 0.25
subplot_bounds3 = [0.05,0.35,0.65]
subplot_diff4 = 0.22
subplot_bounds4 = [0.05,0.275,0.5,0.725]

F1 = aplpy.FITSFigure(inth2co.hdu, subplot=[0.05,0.1,0.25,0.8], figure=fig)
F2 = aplpy.FITSFigure(pcainth2co.hdu, subplot=[0.35,0.1,0.25,0.8], figure=fig)
F3 = aplpy.FITSFigure(pca1inth2co.hdu, subplot=[0.65,0.1,0.25,0.8], figure=fig)


for F in (F1,F2,F3):
    F.show_grayscale()
    F.add_colorbar()
    F.set_tick_labels_format('d.dd','d.dd')
    F.set_tick_xspacing(0.1)
    
for F in (F2,F3,):
    F.tick_labels.hide_y()
    F.axis_labels.hide_y()
F1._ax1.set_title("Before PCA")
F2._ax1.set_title("After PCA-3")
F3._ax1.set_title("After PCA-1")


fig.savefig("Before_vs_After_PCA_clean_3panel.pdf", bbox_inches='tight', dpi=144)


fig = pl.figure(figsize=(16,8))

F1 = aplpy.FITSFigure(diffh2co, subplot=[0.05,0.1,0.4,0.8], figure=fig)
F2 = aplpy.FITSFigure(diffh2co2, subplot=[0.5,0.1,0.4,0.8], figure=fig)

for F in (F1,F2,):
    F.show_grayscale()
    F.add_colorbar()
    F.set_tick_labels_format('d.dd','d.dd')
    F.set_tick_xspacing(0.1)
    
for F in (F2,):
    F.tick_labels.hide_y()
    F.axis_labels.hide_y()
F1._ax1.set_title("raw-PCA3")
F2._ax1.set_title("raw-PCA1")


fig.savefig("PCA_Diffs_3pca_1pca_2panel.pdf", bbox_inches='tight', dpi=144)

spec = h2cocube.mean(axis=(1,2))
pcaspec = pcacube.mean(axis=(1,2))
pca1spec = pca1cube.mean(axis=(1,2))



pl.figure()
pl.subplot(2,2,1)
pl.plot(spec, linewidth=0.5)
pl.plot((spec-pcaspec))
pl.ylim(-0.005,0.007)
pl.xlim(1500,2500)

pl.subplot(2,2,3)
pl.plot(spec, linewidth=0.5)
pl.plot((spec-pca1spec))
pl.ylim(-0.005,0.007)
pl.xlim(1500,2500)


#pl.subplot(3,2,5)
#pl.plot(spec.value, linewidth=0.5)
#pl.plot((spec-timepcaspec).value)
#pl.ylim(-0.005,0.007)
#pl.xlim(1500,2500)

pl.subplot(2,2,2)
pl.plot(spec, linewidth=0.5)
pl.plot((pcaspec), alpha=0.5)
pl.ylim(-0.005,0.007)
pl.xlim(1500,2500)

pl.subplot(2,2,4)
pl.plot(spec, linewidth=0.5)
pl.plot((pca1spec), alpha=0.5)
pl.ylim(-0.005,0.007)
pl.xlim(1500,2500)


