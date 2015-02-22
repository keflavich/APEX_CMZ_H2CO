Data Table README
=================

Column Descriptions
-------------------

fitted_line_parameters_Chi2Constraints.ipac:

Source_Name                 str       Source name in the accompanying .reg files
ComponentID                 str       Component # for the source (there may be multiple velocity components per source)
GLON                        str       Galactic longitude
GLAT                        str       Galactic latitude
ampH2CO                     str       Fitted H2CO 303 amplitude (T_A*, K)
eampH2CO                    str       Error on ampH2CO
h2coratio321303             str       Fitted 303/321 ratio
eh2coratio321303            str       Error on h2coratio321303
center                      str       Velocity centroid (LSR, km/s)
ecenter                     str       Error on center
width                       str       Velocity width (sigma, km/s)
ewidth                      str       Error on width
h2coratio322321             str       Fitted 321/322 ratio
eh2coratio322321            str       Error on h2coratio322321
ampCH3OH                    str       Fitted CH3OH amplitude (T_A*, K)
eampCH3OH                   str       Error on ampCH3OH
spline_ampCH3OH             str       Spline-baselined spectra: Fitted H2CO 303 amplitude (T_A*, K)
espline_ampCH3OH            str       Spline-baselined spectra: Error on ampH2CO
spline_center               str       Spline-baselined spectra: Fitted 303/321 ratio
espline_center              str       Spline-baselined spectra: Error on h2coratio321303
spline_h2coratio321303      str       Spline-baselined spectra: Velocity centroid (LSR, km/s)
espline_h2coratio321303     str       Spline-baselined spectra: Error on center
spline_width                str       Spline-baselined spectra: Velocity width (sigma, km/s)
espline_width               str       Spline-baselined spectra: Error on width
spline_h2coratio322321      str       Spline-baselined spectra: Fitted 321/322 ratio
espline_h2coratio322321     str       Spline-baselined spectra: Error on h2coratio322321
spline_ampH2CO              str       Spline-baselined spectra: Fitted CH3OH amplitude (T_A*, K)
espline_ampH2CO             str       Spline-baselined spectra: Error on ampCH3OH
boxwidth                    str       Width of the box region (deg)
boxheight                   str       Height of the box region (deg)
radius                      str       Radius of the circular region (deg)
area                        str       Region area (deg^2)
posang                      str       Position angle of the box region
is_good                     str       Flag: Are the line fits good enough to use in plotting?
higalcolumndens             str       Herschel Hi-Gal derived mean column density (cm^-2)
higaldusttem                str       Herschel Hi-Gal derived mean temperature (K)
temperature_chi2            str       Chi^2 fitting derived gas temperature (K)
tmin1sig_chi2               str       Chi^2 1-sigma lower limit on temperature 
tmax1sig_chi2               str       Chi^2 1-sigma upper limit on temperature 
column_chi2                 str       Chi^2 fitting derived H2CO column (cm^-2 / (km/s/pc))
cmin1sig_chi2               str       Chi^2 1-sigma lower limit on column
cmax1sig_chi2               str       Chi^2 1-sigma upper limit on column
density_chi2                str       Chi^2 fitting derived gas density (cm^-3)
dmin1sig_chi2               str       Chi^2 1-sigma lower limit on density
dmax1sig_chi2               str       Chi^2 1-sigma upper limit on density
logh2column                 str       Log of the Herschel hi-gal column density (log[cm^-2])
elogh2column                str       
logabundance                str       Log abundance of H2CO relative to H2
elogabundance               str        
tkin_turb                   str       Predicted kinetic temperature of the gas from turbulent heating (K)
reff_pc                     str       Effective radius of the region (pc)
