#!/usr/local/bin/python

'''
MOSFIRE 'handle' command:

(c) npk - Dec 2013
'''
import MOSFIRE
import MOSFIRE.IO as IO
import os
import numpy as np
import pyfits as pf
import sys
import glob



if len(sys.argv) < 3:
    print '''Usage: mospy handle [target]'''
    sys.exit()

files = []
for i in range(1, len(sys.argv)):
    files.extend(glob.iglob(sys.argv[i]))

masks = {}

for fname in files:

    try:
        header = MOSFIRE.IO.readheader(fname)
    except IOError, err:
        print "Couldn't IO %s" % fname
        continue
    except:
        print "%s is unreadable" % fname
        continue

    lamps = ""
    try:
        if header["pwstata7"] == 1:
            lamps += header["pwloca7"][0:2]
        if header["pwstata8"] == 1:
            lamps += header["pwloca8"][0:2]
    except KeyError:
        lamps = "???"
        
    header.update("lamps", lamps)

    try:
        if header["aborted"]:
            header.update("object", "ABORTED")
    except:
        print "Missing header file in: %s" % fname

    try:
        print "%(datafile)12s %(object)40s %(truitime)6.1f s %(maskname)35s %(lamps)3s %(filter)4s %(mgtname)7s" % (header)
    except:
        try:
            print "%(datafile)12s %(object)25s %(truitime)6.1f s %(lamps)3s %(filter)6s %(mgtname)7s" % (header)
        except:
            print "%s Skipped" % fname
            continue


    datafile = header['datafile'] + '.fits'
    maskname = header['maskname']
    filter = header['filter']
    yr,mn,dy = IO.fname_to_date_tuple(datafile)
    date = str(yr)+mn+str(dy)
    object = header['object']

    itime = header['truitime']
    grating_turret = header['mgtname']

    if object.find("MIRA") == -1: mira = False
    else: mira = True

    if maskname.find(" (align)") == -1:
        align = False
    else:
        maskname = maskname.replace(" (align)", "")
        align = True

    if maskname.find('LONGSLIT') != -1:
        align = False

    if maskname.find('long2pos') != -1:
        if grating_turret != 'mirror':
            align = False

    empty_files = {'Align': [], 'Ne': [], 'Ar': [], 'Flat': [],
            'Dark': [], 'Aborted': [], 'Image': [], 'MIRA': [], 'Unknown': []}

    if maskname not in masks:
        masks[maskname] = {date: {filter: empty_files}}

    if date not in masks[maskname]:
        masks[maskname][date] = {filter: empty_files}

    if filter not in masks[maskname][date]:
        masks[maskname][date][filter] = empty_files


    offset = 'Offset_' + str(header['YOFFSET'])

    

    if mira:
        masks[maskname][date][filter]['MIRA'].append(fname)
    elif align:
        masks[maskname][date][filter]['Align'].append(fname)
    elif 'Ne' in header['lamps']:
        masks[maskname][date][filter]['Ne'].append(fname)
    elif 'Ar' in header['lamps']:
        masks[maskname][date][filter]['Ar'].append(fname)
    elif header['ABORTED']:
        masks[maskname][date][filter]['Aborted'].append(fname)
    elif header['FILTER'] == 'Dark':
        masks[maskname][date][filter]['Dark'].append(fname)
    elif header['FLATSPEC'] == 1:
        masks[maskname][date][filter]['Flat'].append(fname)
    elif header['mgtname'] == 'mirror':
        masks[maskname][date][filter]['Image'].append(fname)
    elif offset != 0:
        if offset in masks[maskname][date][filter]: 
            masks[maskname][date][filter][offset].append((fname, itime))
        else: 
            masks[maskname][date][filter][offset] = [(fname, itime)]
    else:
        masks[maskname][date][filter]['Unknown'].append(fname)



##### Now handle mask dictionary

def descriptive_blurb():
    import getpass, time

    uid = getpass.getuser()
    date = time.asctime()

    return "# Created by '%s' on %s\n" % (uid, date)


# Write out the list of files in filepath
#   list = ['/path/to/mYYmmDD_####.fits' ...]
#   filepath is absolute path to the file name to write to
#
#   Result, is a file called filepath is written with
# fits files in the list.
def handle_file_list(output_file, files):
    '''Write a list of paths to MOSFIRE file to output_file.'''

    if os.path.isfile(output_file):
        print "%s: already exists, skipping" % output_file 
        pass

    print "\t", output_file
    f = open(output_file, "w")
    f.write(descriptive_blurb())
    if len(files) == 0:
        f.close()
        return

    picker = lambda x: x
    if len(files[0]) == 2: picker = lambda x: x[0]

    # Identify unique path to files:
    paths = [os.path.dirname(picker(file)) for file in files]
    paths = list(set(paths))

    if len(paths) == 1:
        path_to_all = paths[0]
        converter = os.path.basename
        f.write("%s # Abs. path to files [optional]\n" % path_to_all)
    else:
        converter = lambda x: x


    for path in files:
        if len(path) == 2:  to_write = "%s # %s s\n" % (converter(path[0]), path[1])
        else:               to_write = "%s\n" % converter(path)

        f.write("%s" % to_write)
            

    f.close()

def handle_date_and_filter(mask, date, filter, mask_info):

    path = os.path.join(mask,date,filter)
    try: os.makedirs(path)
    except OSError: pass

    for type in mask_info.keys():
        handle_file_list(os.path.join(path, type + ".txt"), mask_info[type])


for mask in masks.keys():
    for date in masks[mask].keys():
        for filter in masks[mask][date].keys():
            handle_date_and_filter(mask, date, filter, masks[mask][date][filter])


