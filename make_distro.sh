#!/bin/csh


set DT="`date +'%Y.%m.%d'`"
echo $DT

rm repository_version
hg sum > repository_version
/usr/local/EPD/epd-7.2.2/bin/python update_options.py $DT
mkdir badpixels
cp /scr2/mosfire/badpixels/badpix_10sep2012.fits badpixels
tar -cfmosdrp_$DT.tar MOSFIRE/*py apps/*py apps/mospy platescale/10March2011.4.972.db badpixels/badpix_10sep2012.fits repository_version drivers/*py

gzip mosdrp_$DT.tar
cp mosdrp_$DT.tar.gz ~/public_html/mosdrp_releases

set PD=`pwd`
cd ~/public_html/mosdrp_releases
rm mosdrp_cur.tgz
ln -s mosdrp_$DT.tar.gz mosdrp_cur.tgz

cd $PD



