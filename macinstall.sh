#!/bin/csh -f

if(`whoami` == 'root') then
    echo "Please execute this as your user account and not in sudo"
    exit
endif

source /usr/stsci/envconfig.mac/cshrc
iraf

/usr/stsci/pyssg/Python-2.7.3/bin/python install.py

rehash

