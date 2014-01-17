#!/bin/csh


set DT="`date +'%Y.%m.%d_%H.%M.%S'`"
echo $DT

mkdir badpixels
cp /scr2/mosfire/badpixels/badpix_10sep2012.fits badpixels
tar -cfmosfire_v$DT.tar MOSFIRE/*py apps/*py apps/mospy platescale/10March2011.4.972.db badpixels/badpix_10sep2012.fits

gzip mosfire_v$DT.tar
