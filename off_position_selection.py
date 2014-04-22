"""
Quick-look to demonstrate that the off positions are relatively emission-free
"""
import aplpy
from astropy.utils.data import download_file
dame2001 = download_file('http://www.cfa.harvard.edu/mmw/Wco_DHT2001.fits.gz', cache=True)
F = aplpy.FITSFigure(dame2001, convention='calabretta')
F.show_grayscale()
F.recenter(0,0,width=5,height=3) # center on the inner 5x3 degrees
F.show_regions('target_fields_8x8.reg')
F.show_regions('off_positions_selectedfromDame2001.reg')
F.save('Dame2001_APEXCMZ_offpositions.png')
F.save('Dame2001_APEXCMZ_offpositions.pdf')
