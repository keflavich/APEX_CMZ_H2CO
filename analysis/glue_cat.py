from glue.core import Data, DataCollection
from glue.core.data_factories import astropy_tabular_data, load_data
from glue.core.link_helpers import LinkSame, LinkTwoWay
from glue.qt.glue_application import GlueApplication
from glue.qt.widgets import ScatterWidget, ImageWidget
from glue.qt.widgets.image_widget import (_slice_from_path, _slice_label,
                                          _slice_index, StandaloneImageWidget,
                                          PVSliceWidget)
from glue import qglue

from astropy import units as u
from astropy import coordinates
from astropy import wcs
from astropy.io import ascii
from paths import mpath,apath,fpath,molpath,hpath
import pvextractor


catalog = astropy_tabular_data(apath('fitted_line_parameters_Chi2Constraints.ipac'),
                               format='ascii.ipac')
catalog.label='FittedLineParameters'
catalog.style.color = 'r'
catalog.style.marker = '+'
cube = load_data('/Users/adam/work/gc/apex/h2co_cubes/APEX_H2CO_303_202_bl.fits')

dc = DataCollection([cube, catalog])

dc.add_link(LinkSame(cube.id['Galactic Longitude'], catalog.id['GLON']))
dc.add_link(LinkSame(cube.id['Galactic Latitude'], catalog.id['GLAT']))

def ms_to_kms(x): return x/1e3
def kms_to_ms(x): return x*1e3

dc.add_link(LinkTwoWay(cube.id['Vrad'], catalog.id['center'], ms_to_kms, kms_to_ms))

subset_tem_lt_60 = (catalog.id['temperature_chi2'] < 60) & (catalog.id['temperature_chi2'] > 10)

subset_tem_gt_60 = catalog.id['temperature_chi2'] > 60

dc.new_subset_group(label='T < 60', subset_state=subset_tem_lt_60)
dc.new_subset_group(label='T > 60', subset_state=subset_tem_gt_60)

dc.subset_groups[0].style.markersize=15
dc.subset_groups[0].style.marker='+'
dc.subset_groups[1].style.markersize=15
dc.subset_groups[1].style.marker='*'

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
coords = coordinates.SkyCoord(table['l']*u.deg, table['b']*u.deg, frame='galactic')[:10]
P = pvextractor.Path(coords)
pv = pvextractor.extract_pv_slice(cube.data['PRIMARY'], P, wcs=cube.data.coords.wcs)
pvwidget = PVSliceWidget(image=pv.data, wcs=wcs.WCS(pv.header),
                         image_widget=cube_viewer, interpolation='nearest')


app.start()
