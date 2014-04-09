import sys


ver = sys.argv[1]

if len (ver.split('.')) != 3:
    raise Exception("Version doesn't seem right")


f = open("MOSFIRE/Options.py")
lines = f.readlines()
f.close()

found = False
for i in xrange(len(lines)):
    line = lines[i]
    sp = line.split("=")
    if len(sp) == 2:
        if sp[0].rstrip() == '__version__':
            sp[1] = " '" + sys.argv[1] + "'\n"

        lines[i] = "=".join(sp)
        found = True
        break

if not found:
    raise Exception("Version not updated")

f = open("MOSFIRE/Options.py", "w")
f.writelines(lines)
f.close()



