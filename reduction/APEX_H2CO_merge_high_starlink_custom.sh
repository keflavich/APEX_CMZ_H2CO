#!/bin/bash
. /star/etc/profile
kappa > /dev/null
convert > /dev/null
fits2ndf APEX_H2CO_merge_high.fits APEX_H2CO_merge_high.sdf
mfittrend APEX_H2CO_merge_high.sdf  ranges=\"218000000000.00 2.181266E+11 2.183057E+11 2.183697E+11 2.184891E+11 2.186718E+11 2.187722E+11 2.188428E+11 2.189484E+11 218999320002.91\" order=2 axis=3 out=APEX_H2CO_merge_high_baseline.sdf
sub APEX_H2CO_merge_high.sdf APEX_H2CO_merge_high_baseline.sdf APEX_H2CO_merge_high_sub.sdf
sqorst APEX_H2CO_merge_high_sub.sdf mode=factors  axis=3 factors=0.3 out=APEX_H2CO_merge_high_vrebin
gausmooth APEX_H2CO_merge_high_vrebin fwhm=2.0 axes=[1,2] out=APEX_H2CO_merge_high_smooth
#collapse APEX_H2CO_merge_high estimator=mean axis="RADI-LSR" low=-400 high=500 out=APEX_H2CO_merge_high_continuum
rm APEX_H2CO_merge_high_sub.fits
ndf2fits APEX_H2CO_merge_high_sub APEX_H2CO_merge_high_sub.fits
rm APEX_H2CO_merge_high_smooth.fits
ndf2fits APEX_H2CO_merge_high_smooth APEX_H2CO_merge_high_smooth.fits
# Fix STARLINK's failure to respect header keywords.
sethead APEX_H2CO_merge_high_smooth.fits RESTFRQ=`gethead RESTFRQ APEX_H2CO_merge_high.fits`

mfittrend APEX_H2CO_merge_high.sdf  ranges=\"218000000000.00 2.181266E+11 2.183057E+11 2.183697E+11 2.184891E+11 2.186718E+11 2.187722E+11 2.188428E+11 2.189484E+11 218999320002.91\" order=5 axis=3 out=APEX_H2CO_merge_high_baseline2.sdf
sub APEX_H2CO_merge_high.sdf APEX_H2CO_merge_high_baseline2.sdf APEX_H2CO_merge_high_sub2.sdf
sqorst APEX_H2CO_merge_high_sub2.sdf mode=factors  axis=3 factors=0.3 out=APEX_H2CO_merge_high_vrebin2
gausmooth APEX_H2CO_merge_high_vrebin2 fwhm=5.0 axes=[1,2] out=APEX_H2CO_merge_high_vsmooth
#regrid APEX_H2CO_merge_high_vsmooth APEX_H2CO_merge_high_vsmoothds mapping=! scale=0.2 method=blockave params=2
rm temp.sdf
sqorst APEX_H2CO_merge_high_vsmooth mode=factors  axis=1 factors=0.4 out=temp
rm APEX_H2CO_merge_high_vsmoothds.sdf
sqorst temp mode=factors  axis=2 factors=0.4 out=APEX_H2CO_merge_high_vsmoothds
rm APEX_H2CO_merge_high_vsmoothds.fits
ndf2fits APEX_H2CO_merge_high_vsmoothds APEX_H2CO_merge_high_vsmoothds.fits

rm APEX_H2CO_merge_high.sdf
rm APEX_H2CO_merge_high_smooth.sdf
rm APEX_H2CO_merge_high_vrebin.sdf
rm APEX_H2CO_merge_high_vsmoothds.sdf
rm APEX_H2CO_merge_high_baseline.sdf
rm APEX_H2CO_merge_high_sub.sdf
rm APEX_H2CO_merge_high_vrebin2.sdf  temp.sdf
rm APEX_H2CO_merge_high_baseline2.sdf
rm APEX_H2CO_merge_high_sub2.sdf
rm APEX_H2CO_merge_high_vsmooth.sdf
