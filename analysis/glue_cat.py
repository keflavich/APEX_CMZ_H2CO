from glue.core import Data, DataCollection
from glue.core.data_factories import astropy_tabular_data, load_data
from glue.core.link_helpers import LinkSame, LinkTwoWay
from glue.qt.glue_application import GlueApplication
from glue.qt.widgets import ScatterWidget, ImageWidget
from glue.qt.widgets.image_widget import (_slice_from_path, _slice_label,
                                          _slice_index, StandaloneImageWidget,
                                          PVSliceWidget)
from glue import qglue

import matplotlib
import numpy as np
from astropy import units as u
from astropy import coordinates
from astropy import wcs
from astropy.io import ascii
from paths import mpath,apath,fpath,molpath,hpath
import pvextractor

catalog = astropy_tabular_data(apath('fitted_line_parameters_Chi2Constraints.ipac'),
                               format='ascii.ipac')
catalog.label='FittedLineParameters'
catalog.style.color = 'green'
catalog.style.marker = 'o'
cube = load_data(hpath('APEX_H2CO_303_202_bl.fits'))
cube.label='H2CO 303/202'
cube2 = load_data(molpath('APEX_SiO_54.fits'))
cube2.label='SiO'
cube3 = load_data(hpath('APEX_13CO_matched_H2CO.fits'))
cube3.label='13CO'
higaltem = load_data('/Users/adam/work/gc/gcmosaic_temp_conv36.fits')

dc = DataCollection([cube, catalog, cube2, cube3, higaltem])
dc.merge(cube,cube2,cube3)

dc.add_link(LinkSame(cube.id['Galactic Longitude'], catalog.id['GLON']))
dc.add_link(LinkSame(cube.id['Galactic Latitude'], catalog.id['GLAT']))

def ms_to_kms(x): return x/1e3
def kms_to_ms(x): return x*1e3

dc.add_link(LinkTwoWay(cube.id['Vrad'], catalog.id['center'], ms_to_kms, kms_to_ms))

subset_tem_lt_60 = (catalog.id['temperature_chi2'] < 60) & (catalog.id['temperature_chi2'] > 10) & (catalog.id['area'] < 0.015)

subset_tem_gt_60 = (catalog.id['temperature_chi2'] > 60) & (catalog.id['area'] < 0.015)

app = GlueApplication(dc)

# plot x vs y, flip the x axis, log-scale y axis
scatter = app.new_data_viewer(ScatterWidget)
scatter.add_data(catalog)
scatter.yatt = catalog.id['temperature_chi2']
scatter.xatt = catalog.id['higaldusttem']
#scatter.xflip = True
#scatter.ylog = True

cube_viewer = app.new_data_viewer(ImageWidget)
cube_viewer.add_data(cube)
#cube_viewer.add_subset(subset_tem_lt_60)
#cube_viewer.add_subset(subset_tem_gt_60)
#cube_viewer.add_data(catalog)

table = ascii.read(apath('orbit_K14_2.dat'), format='basic', comment="#", guess=False) 
coords = coordinates.SkyCoord(table['l']*u.deg, table['b']*u.deg, frame='galactic')
P = pvextractor.Path(coords, width=120*u.arcsec)

x, y = [np.round(_x).astype(int) for _x in P.sample_points(1, wcs=cube.data.coords.wcs)]

ax = cube_viewer.axes
ax.set_axis_bgcolor('black')
#ax.plot(table['l'], table['b'], 'r-', linewidth=2, alpha=0.5)
#ax.plot(x, y, 'r-', linewidth=2, alpha=0.5)
patches = P.to_patches(1, ec='red', fc='none',
                       #transform=ax.transData,
                       clip_on=True, #clip_box=ax.bbox,
                       wcs=cube.data.coords.wcs)
for patch in patches:
    patch.set_linewidth(0.5)
    patch.set_alpha(0.5)
    patch.zorder = 50
patchcoll = matplotlib.collections.PatchCollection(patches, match_original=True)
patchcoll.zorder=10
ax.add_collection(patchcoll)
ax.axis([x.min(),x.max(),y.min(),y.max()])

pv = pvextractor.extract_pv_slice(cube.data['PRIMARY'], P, wcs=cube.data.coords.wcs)
pvwidget = PVSliceWidget(image=pv.data, x=x, y=y,# wcs=wcs.WCS(pv.header),
                         image_widget=cube_viewer, interpolation='nearest')
pv_viewer = app.add_widget(pvwidget, label="Orbit PV Slice")
ax2 = pvwidget.axes

dl = (table['l'][1:]-table['l'][:-1])
db = (table['b'][1:]-table['b'][:-1])
dist = (dl**2+db**2)**0.5
cdist = np.zeros(dist.size+1) * u.deg
cdist[1:] = dist.cumsum() * u.deg
#pixscale = ((x[1]-x[0])**2+(y[1]-y[0])**2)**0.5
pixscale = wcs.utils.celestial_pixel_scale(cube.data.coords.wcs)
spwcs = cube.data.coords.wcs.sub([wcs.WCSSUB_SPECTRAL])
spax = spwcs.wcs_world2pix(table["v'los"]*1e3, 0)[0]
ax2.plot(cdist/pixscale, spax, 'r-', linewidth=2, alpha=0.5)
ax2.set_axis_bgcolor('black')

dc.new_subset_group(label='T < 60', subset_state=subset_tem_lt_60)
dc.new_subset_group(label='T > 60', subset_state=subset_tem_gt_60)

dc.subset_groups[0].style.markersize=15
dc.subset_groups[0].style.marker='+'
dc.subset_groups[0].style.color='blue'
dc.subset_groups[1].style.markersize=15
dc.subset_groups[1].style.marker='*'
dc.subset_groups[1].style.color='orange'

# SERIOUSLY, DO IT
ax.axis([x.min(),x.max(),y.min(),y.max()])

app.start()
