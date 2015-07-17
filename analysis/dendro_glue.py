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
dendcube.label='Dendrogram Cube' # this label is redundant, it will be deleted upon merge
# cube contains real WCS information; dendcube does not
h2cocube = load_data(hpath('APEX_H2CO_303_202_bl.fits'))
h2cocube.label='H2CO 303202 Cube'
catalog = astropy_tabular_data(hpath('PPV_H2CO_Temperature.ipac'), format='ipac')
catalog.label='Fitted Catalog'


h2cocube.add_component(dendcube.get_component('structure'), 'structure')

dc = DataCollection(dendrogram)
dc.append(h2cocube)
dc.append(catalog)
dc.append(dendcube)

dc.merge(h2cocube,dendcube)
dc.merge(dendro, catalog)

app = GlueApplication(dc)

cube_viewer = app.new_data_viewer(ImageWidget)
cube_viewer.add_data(h2cocube)

h2cocube.join_on_key(dendro, 'structure', dendro.pixel_component_ids[0])

scatter = app.new_data_viewer(ScatterWidget)
scatter.add_data(dendro)
scatter.yatt = dendro.id['temperature_chi2']
scatter.xatt = catalog.id['r321303']

dendview = app.new_data_viewer(DendroWidget)
dendview.add_data(dendro)

subset_tem_bt_40_60 = ((catalog.id['temperature_chi2'] < 60) &
                    (catalog.id['temperature_chi2'] > 40))
subset_tem_lt_40 = ((catalog.id['temperature_chi2'] < 40) &
                    (catalog.id['temperature_chi2'] > 10))
subset_tem_gt_60 = (catalog.id['temperature_chi2'] > 60)
dc.new_subset_group(label='T < 40', subset_state=subset_tem_lt_40)
dc.new_subset_group(label='40 < T < 60', subset_state=subset_tem_bt_40_60)
dc.new_subset_group(label='T > 60', subset_state=subset_tem_gt_60)

#start Glue
app.start()
