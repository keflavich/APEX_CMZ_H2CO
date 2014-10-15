"""
Quick-look to demonstrate that the off positions are relatively emission-free
"""
import aplpy
from astropy.utils.data import download_file
import paths
dame2001 = download_file('http://www.cfa.harvard.edu/mmw/Wco_DHT2001.fits.gz', cache=True)
F = aplpy.FITSFigure(dame2001, convention='calabretta')
F.show_grayscale()
F.recenter(0,0,width=5,height=3) # center on the inner 5x3 degrees
F.show_regions(paths.rpath('target_fields_8x8.reg'))
F.show_regions(paths.rpath('off_positions_selectedfromDame2001.reg'))
F.save(paths.fpath('Dame2001_APEXCMZ_offpositions.png'))
F.save(paths.fpath('Dame2001_APEXCMZ_offpositions.pdf'))

F.hide_layer('region_set_1')
F.hide_layer('region_set_1_txt')
F.show_regions(paths.rpath('target_fields_8x8_coloredbyoffposition.reg'))
F.save(paths.fpath('Dame2001_APEXCMZ_offpositions_coloredbyoff.png'))
