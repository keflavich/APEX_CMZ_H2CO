#!/bin/csh -f
#
#

set newname = $1    
echo $newname
if ($newname == "") then
    echo "Usage: ./rename.csh <PROJECT_ID>"
    echo "e.g. ./rename.csh M0023_93"
    exit
endif

set base = `echo $newname | awk -F"_" '{print $1}' | sed "s/M//"`
echo $base
set period = `echo $newname | awk -F"_" '{print $2}'`
echo $period

@ year = ( $period + 1 ) / 2 + 1967
echo $year

set newid = "M-0"$period".F-"$base"-"$year
echo $newid

foreach file ( M0000_00* )
  set newfile = `echo $file | sed "s/M0000_00/$newname/g"`
  sed "s/M0000_00/$newname/g" $file | sed "s/M-000.F-0000-0000/$newid/g" > $newfile
end

