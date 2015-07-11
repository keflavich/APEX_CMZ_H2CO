"""
Load dendrograms into a Glue session
"""

from glue.core.data_factories import load_data
try:
    # v>=0.5
    from glue.core.data_factories.dendro_loader import load_dendro
    from glue.core.data_factories.tables import astropy_tabular_data
except ImportError:
    # v<=0.4
    from glue.core.data_factories import load_dendro
    from glue.core.data_factories import astropy_tabular_data
from glue.core import DataCollection
from glue.core.link_helpers import LinkSame
from glue.qt.glue_application import GlueApplication
from glue.core import Data, DataCollection, Component
from glue.core.link_helpers import LinkSame, LinkTwoWay
from glue.qt.glue_application import GlueApplication
from glue.qt.widgets import ScatterWidget, ImageWidget
from glue.qt.widgets.dendro_widget import DendroWidget
from glue.qt.widgets.image_widget import StandaloneImageWidget
from glue import qglue

import matplotlib
import numpy as np
from astropy import log
from astropy import units as u
from astropy import coordinates
from astropy import wcs
from astropy.table import Table
from astropy.io import ascii
try:
    from paths import mpath,apath,fpath,molpath,hpath
except ImportError:
    hpath = lambda x:x


dendrogram = load_dendro(hpath('DendroMask_H2CO303202.hdf5'))
dendro,dendcube = dendrogram
dendcube.label='Dendrogram Cube'
# cube contains real WCS information; dendcube does not
h2cocube = load_data(hpath('APEX_H2CO_303_202_bl.fits'))
h2cocube.label='H2CO 303202 Cube'
catalog = astropy_tabular_data(hpath('PPV_H2CO_Temperature.ipac'), format='ipac')
catalog.label='Fitted Catalog'


print()
log.info("cube components: {}".format(h2cocube.components))
log.info("dendcube components: {}".format(dendcube.components))
log.info("catalog components: {}".format(catalog.components))
log.info("dendro components: {}".format(dendro.components))

h2cocube.add_component(dendcube.get_component('structure'), 'structure')
h2cocube.join_on_key(dendro, 'structure', dendro.pixel_component_ids[0])

dc = DataCollection([h2cocube, dendrogram, catalog, dendcube])

# not obvious whether any merging is needed
#dc.merge(h2cocube,dendcube)
#dc.merge(catalog, dendro)

dc.append(catalog)

app = GlueApplication(dc)

cube_viewer = app.new_data_viewer(ImageWidget)
cube_viewer.add_data(h2cocube)

# not obvious whether these are necessary:
dc.add_link(LinkSame(h2cocube.id['Pixel x'], dendcube.id['Pixel x']))
dc.add_link(LinkSame(h2cocube.id['Pixel y'], dendcube.id['Pixel y']))
dc.add_link(LinkSame(h2cocube.id['Pixel z'], dendcube.id['Pixel z']))

dc.add_link(LinkSame(h2cocube.id['Galactic Longitude'], dendcube.id['Galactic Longitude']))
dc.add_link(LinkSame(h2cocube.id['Galactic Latitude'], dendcube.id['Galactic Latitude']))
dc.add_link(LinkSame(h2cocube.id['Vrad'], dendcube.id['Vrad']))

dc.add_link(LinkSame(catalog.id['_idx'], dendro.id['index']))

def ms_to_kms(x): return x/1e3
def kms_to_ms(x): return x*1e3

scatter = app.new_data_viewer(ScatterWidget)
scatter.add_data(catalog)
scatter.yatt = catalog.id['temperature_chi2']
#scatter.xatt = catalog.id['r321303']
scatter.xatt = catalog.id['_idx']

# selection on the '_idx' parameter works, but absolutely nothing else does.

dendview = app.new_data_viewer(DendroWidget)
dendview.add_data(dendro)

#start Glue
app.start()
