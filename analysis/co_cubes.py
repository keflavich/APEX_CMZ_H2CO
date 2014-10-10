from spectral_cube import SpectralCube
from paths import hpath

cube13co = SpectralCube.read(hpath('APEX_13CO_matched_H2CO.fits'))
cube18co = SpectralCube.read(hpath('APEX_C18O_matched_H2CO.fits'))
cube13cosm = SpectralCube.read(hpath('APEX_13CO_matched_H2CO_smooth.fits'))
cube18cosm = SpectralCube.read(hpath('APEX_C18O_matched_H2CO_smooth.fits'))
