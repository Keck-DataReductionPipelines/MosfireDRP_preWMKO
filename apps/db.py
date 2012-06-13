

import os
import pdb
import sqlite3
import sys

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

        if len(path.split("/")) > 1: continue

        try: date = int(path[0:4])
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

commands = [make]
def usage():
    print """
Commands: """

    for command in commands:
        print("\t" + command.__name__ + ": " + command.__doc__)

    print("\n")

if __name__ == '__main__':

    if len(sys.argv) != 3:
        usage()
        sys.exit()

    if sys.argv[2] == 'make':
        print "Making database"
        make()
    else:
        usage()
        sys.exit()
