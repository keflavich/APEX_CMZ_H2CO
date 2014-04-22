APEX_CMZ_H2CO
=============

Observing scripts and logs from the 2013-2014 APEX survey of 218 GHz H2CO in
the CMZ 

The MPI and ESO observing scripts are in their respective directories.  A few
"emission-free" off positions were selected based on the Dame 2001 data.
Run [off_position_selection.py](off_position_selection.py) to generate a diagram
demonstrating this (it will automatically download the necessary data).  One
can also demonstrate that these are relatively emission-free using the [Bell
Labs 7m 12CO 1-0 data](http://files.figshare.com/1216354/GC_12CO_LVcube.fits),
with apologies for the FITS header.

The off position map, with two additional points identified by Jens Kauffmann
and Walker Lu:
![Imgur](http://i.imgur.com/Oh1HI1v.png)

Note that Jens' point reportedly showed some signal around SiO, while Walker's
did not.
