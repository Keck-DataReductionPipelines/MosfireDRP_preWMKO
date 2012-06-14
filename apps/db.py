

import os
import pdb
import sqlite3
import sys

from operator import itemgetter
from itertools import groupby

import MOSFIRE 

from MOSFIRE import Options, IO




def load_db():
    indir = Options.indir
    outname = os.path.join(Options.outdir, "mosfire_files.db")

    print("Database located at {0}".format(outname))
    conn = sqlite3.connect(outname)

    return conn

def create(cursor):
    cursor.execute('''
        CREATE TABLE if not exists files
        (id integer primary key, path text, fdate text, number integer)
    ''')

keys = []

def append_column(cursor, name, typename):
    qry = "alter table files\nadd {0} {1}".format(name, typename)
    try:
        cursor.execute(qry)
        print("Added {0} as {1}".format(name, typename))
    except sqlite3.OperationalError:
        pass


def make():
    """Make the database"""

    db = load_db()
    c = db.cursor()
    create(c)
    dirs = os.walk(Options.indir)

    for root, dirs, files in dirs:
        if root == Options.indir: continue
        ignore, path = root.split(Options.indir)

        if len(path.split("/")) != 2: continue

        try: date = int(path.split("/")[1][0:4])
        except: continue

        if (date < 2012) or (date > 2030): continue



        for file in files:
            if len(file) != 17: continue
            p = os.path.join(root, file)

            num = db.execute('select count(*) from files where path = "%s"' %
                    p).fetchall()
            if num[0][0] > 0: 
                print("Skipping: " + p)
                continue
            print(p)

            hdr = IO.readheader(p)

            fdate = file.split("_")[0][1:]
            number = file.split("_")[1][:-5]
            insert_sql = "insert into files(path,fdate,number,"
            vals = "?,?,?,"
            values = [p, fdate, number]

            for key in hdr.keys():
                value = hdr[key]
                T = type(value)
                key = key.replace("-","_")

                insert_sql += key + ","
                vals += "?,"
                values.append(value)

                if key in keys: continue
                keys.append(key)

                if T == int: typename = 'integer'
                if T == float: typename = 'real'
                else: typename = 'text'
                append_column(c, key, typename)


            insert_sql = insert_sql[:-1] + ") values (" + vals[:-1] + ")"
            c.execute(insert_sql, tuple(values))

            

    db.commit()

def find_continuous(data):
    '''Find all continuous numbers in a list'''
    # http://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
    ranges = []
    for k, g in groupby(enumerate(data), lambda (i,x):i-x):
        group = map(itemgetter(1), g)
        ranges.append((group[0], group[-1]))
    return ranges


      

def masks():
    """List all slit masks"""


    db = load_db()


    if len(sys.argv) == 3:
        cur = db.execute("select maskname, count(maskname) from files group by maskname")
        ress = cur.fetchall()
        print("{0:74s} {1:5s}".format("Mask Name", "Count"))
        print("-"*80)
        for res in ress:
            print("{0:74s} {1:5g}".format(res[0], res[1]))

    if len(sys.argv) == 4:
        maskname = sys.argv[3]
        cur = db.execute(
    '''
    select count(filter), filter, itime/1000.0, yoffset 
    from files 

    where maskname = "{0}" and substr(obsmode, -12, 12) = "spectroscopy"

    group by filter'''.
                format(maskname))

        filters = cur.fetchall()

        for filter in filters:
            print
            print("{0:45s} {1:4s}".format(maskname, filter[1]))
            print("-" * 55)

            cur = db.execute('''
    select path, fdate, number
    from files
    where maskname = "{0}" and substr(obsmode, -12, 12) = "spectroscopy" and
    filter = "{1}" and (itime/1000.0) < 30 and (el-45) < .1
    order by fdate, number
            '''.format(maskname, filter[1]))

            FL = cur.fetchall() 

            print "%i flats on %i nights " % (len(FL), len(set([str(S[1]) for
                S in FL])))


            cur = db.execute('''
    select fdate
    from files
    where maskname = "{0}" and substr(obsmode, -12, 12) = "spectroscopy" and
    filter = "{1}" and (itime/1000.0) > 30
    group by fdate
            '''.format(maskname, filter[1]))

            DATES = cur.fetchall()

            for date in DATES:
                date = date[0]
                cur = db.execute('''
    select path, fdate, number, yoffset, itime/1000.0
    from files
    where maskname = "{0}" and substr(obsmode, -12, 12) = "spectroscopy" and
    filter = "{1}" and (itime/1000.0) > 30 and fdate = {2}
    order by fdate, number
    '''.format(maskname, filter[1], date))
                FRAMES = cur.fetchall()

                print("{0}: {1} frames:".format(date, len(FRAMES)))

                nums = [int(S[2]) for S in FRAMES]
                observations = find_continuous(nums)

                for observation in observations:
                    offsets = {}
                    for frame in FRAMES:
                        path, fdate, number, yoffset, itime = frame
                        if offsets.has_key(yoffset):
                            offsets[yoffset]["paths"].append(path)
                            offsets[yoffset]["itime"] += itime
                        else:
                            offsets[yoffset] = {}
                            offsets[yoffset]["paths"] = [path]
                            offsets[yoffset]["itime"] = itime


                for k,v in offsets.iteritems():
                    print("\tOffset {0:5s} has {1:3g} frames total exptime is "
                        "{2:5g} s".format(k, len(v["paths"]), v["itime"]))

            




                

commands = [make, masks]
def usage():
    print """
Commands: """

    for command in commands:
        print("\t" + command.__name__ + ": " + command.__doc__)

    print("\n")

if __name__ == '__main__':

    if len(sys.argv) < 3:
        usage()
        sys.exit()

    if sys.argv[2] == 'make':
        print "Making database"
        make()
    if sys.argv[2] == 'masks':
        masks()

    else:
        usage()
        sys.exit()
