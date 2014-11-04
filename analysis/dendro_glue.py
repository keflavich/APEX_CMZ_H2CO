"""
Load dendrograms into a Glue session
"""

from glue.core.data_factories import load_data,load_dendro
from glue.core import DataCollection
from glue.core.link_helpers import LinkSame
from glue.qt.glue_application import GlueApplication
from glue.core import Data, DataCollection, Component
from glue.core.data_factories import astropy_tabular_data, load_data
from glue.core.link_helpers import LinkSame, LinkTwoWay
from glue.qt.glue_application import GlueApplication
from glue.qt.widgets import ScatterWidget, ImageWidget
from glue.qt.widgets.dendro_widget import DendroWidget
from glue.qt.widgets.image_widget import StandaloneImageWidget
from glue import qglue

import matplotlib
import numpy as np
from astropy import units as u
from astropy import coordinates
from astropy import wcs
from astropy.table import Table
from astropy.io import ascii
try:
    from paths import mpath,apath,fpath,molpath,hpath
except ImportError:
    hpath = lambda x:x


#load 2 datasets from files
dendrogram = load_dendro(hpath('DendroMask_H2CO303202.hdf5'))
dendro,sncube = dendrogram
sncube.label='S/N Cube'
cube = load_data(hpath('APEX_H2CO_303_202_bl.fits'))
table = ascii.read(hpath('PPV_H2CO_Temperature.ipac'), format='ipac')
table['glon'] = table['lon'] - 360*(table['lon'] > 180)
table['xpix'] = table['x_cen'] # Glue "eats" these
table['ypix'] = table['y_cen'] # Glue "eats" these

catalog=Data(parent=table['parent'], label='Fitted Catalog')
#catalog=Data()
for column_name in table.columns:
    cc = table[column_name]
    uu = cc.unit if hasattr(cc, 'unit') else cc.units
    if cc.name == 'parent':
        cc.name = 'cat_parent'
        column_name = 'cat_parent'
    elif cc.name == 'height':
        cc.name = 'cat_height'
        column_name = 'cat_height'
    elif cc.name == 'peak':
        cc.name = 'cat_peak'
        column_name = 'cat_peak'

    nc = Component.autotyped(cc, units=uu)
    catalog.add_component(nc, column_name)
    #  if column_name != 'parent' else '_flarent_'


catalog.join_on_key(dendro, '_idx', dendro.pixel_component_ids[0])
dc = DataCollection(dendrogram)
#dc = DataCollection([cube, dendrogram, catalog])
#dc.merge(cube,sncube)
#sncube.join_on_key(dendro, 'structure', dendro.pixel_component_ids[0])
#dc.merge(catalog, dendro)

# UNCOMMENT THIS LINE TO BREAK THE VIEWER
dc.append(catalog)

app = GlueApplication(dc)

cube_viewer = app.new_data_viewer(ImageWidget)
cube_viewer.add_data(sncube)

# link positional information
dc.add_link(LinkSame(sncube.id['structure'], catalog.id['_idx']))
#dc.add_link(LinkSame(image.id['World y: DEC--TAN'], catalog.id['DEJ2000']))

dc.add_link(LinkSame(cube.id['Galactic Longitude'], catalog.id['x_cen']))
dc.add_link(LinkSame(cube.id['Galactic Latitude'], catalog.id['y_cen']))

def ms_to_kms(x): return x/1e3
def kms_to_ms(x): return x*1e3

dc.add_link(LinkTwoWay(cube.id['Vrad'], catalog.id['v_cen'], ms_to_kms, kms_to_ms))

scatter = app.new_data_viewer(ScatterWidget)
scatter.add_data(catalog)
scatter.yatt = catalog.id['temperature_chi2']
scatter.xatt = catalog.id['r303321']

dendview = app.new_data_viewer(DendroWidget)
dendview.add_data(dendro)

#start Glue
app.start()
