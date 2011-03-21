
import numpy as np
import pylab as pl
import scipy as sp
import scipy.io
import os

path = "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" 
proto = "m11031%1.1i_%4.4i.fits.sav.mat"

pl.ion()
rs = []
exclude = [266, 391,392,393,394,426]
#for i in range(232,554):
for i in range(232,554):
        if i in exclude: continue
        if os.path.exists(path + proto % (2, i)):
                ld = scipy.io.loadmat(path + proto % (2, i))
        elif os.path.exists(path + proto % (3, i)):
                ld = scipy.io.loadmat(path + proto % (3, i))
        else:
                print i
                raise Exception("No luck")

        ld["img_num"] = i
        # Hack for b4
        ld["deltas"][3] += 5
        ld["poss"][3] = ld["poss"][5]
        ld["request"][3] += 2

        rs.append(ld)



def round100(n):
        return np.round(n/30.)*30.

def find_unique_pos(res):

        pl.figure(3)
        pl.clf()
        poss = []
        for r in rs:

                top = np.unique(round100(r["poss"]).ravel())
                top = top[np.isfinite(top)]
                poss.extend(top)

                ok = np.isfinite(r["poss"])
                pl.plot(r["request"][ok], r["deltas"][ok], '*')

                bad = np.abs(r["poss"] - 1901) < .5
                if bad.any(): 
                        print "p: ", r["img_num"], r["poss"][bad].ravel(), r["bars"][bad].ravel(), r["request"][bad]

                bad = (r["deltas"] < -5) | (r["deltas"] > 1)
                if bad.any(): 
                        print "D: ", r["img_num"], r["poss"][bad].ravel(), r["bars"][bad].ravel(), r["deltas"][bad].ravel()

                bad = np.where(np.abs(r["request"]-28.8) < .5)[0]
                bad.ravel()
                if bad.any():
                        print "Last: ", r["img_num"], r["poss"][bad].ravel()

                bad = np.where(np.abs(r["request"]-194.8) < .5)[0]
                bad.ravel()
                if bad.any():
                        print "@: ", r["img_num"], r["poss"][bad].ravel(), r["bars"][bad].ravel(), r["request"][bad+2].ravel(),r["request"][bad].ravel(),r["request"][bad-2].ravel()

        poss = np.unique(np.array(poss))
        ok = np.isfinite(poss)

        return poss[ok]

poss = find_unique_pos(rs)
npos = len(poss)

bars = rs[0]["bars"]
nbars = len(rs[0]["bars"])

def results_to_matrix(rs):
        global npos, nbars, bars, poss

        results = [[ [] for i in range(npos)] for j in range(nbars)]
        for r in rs:
                for bar in range(len(bars)):
                        assert(r["bars"][bar][0] == np.array(bar+1))

                        pos = np.where(round100(r["poss"][bar])  == poss)[0]
                        if len(pos) == 0: continue
                        elif len(pos) > 1:
                                assert(False)

                        results[bar][pos].append(r["deltas"][bar])
        return results
        

results = results_to_matrix(rs)



means = np.zeros((nbars, npos))
stds = np.zeros((nbars, npos))
means[:,:] = np.nan

for i in range(nbars):
        for j in range(npos):
                v = results[i][j]
                if len(v) == 0: continue
                means[i][j] = np.mean(v)
                stds[i][j] = np.std(v)


ss = []
for i in range(1, nbars, 2):
        ok = np.isfinite(means[i,:])
        m = means[i,ok]
        s = stds[i,ok]
        ss.append(np.std(m))

print "RMS Deviation from predicted: %3.2f" % np.mean(ss)
pl.figure(2)
pl.clf()

def draw_bar(xs, y, dat):
        global scale

        for ix in range(len(xs)):
                if dat[ix] > 0: color = "red"
                else: color = "blue"
                pl.arrow(xs[ix], y/2., scale*dat[ix], 0, color=color)


all_shifts = []
def draw_slit(odd, even, iodd):
        ok = np.isfinite(odd)
        shift = odd[ok].mean()
        all_shifts.append(shift)

        odd -= shift
        xodd = poss[ok]
        draw_bar(xodd, iodd, odd[ok])

        ok = np.isfinite(even)
        shift = even[ok].mean()
        all_shifts.append(shift)
        even -= shift
        xeven = poss[ok]
        draw_bar(xeven, iodd+.25, even[ok])
                        
        ok = np.isfinite(even)
        xeven = poss[ok]

scale = 200
for i in range(1, nbars-1, 2):

        odd_deltas = means[i,1:-1]
        even_deltas= means[i+1,2:-1]

        draw_slit(odd_deltas, even_deltas, i)

        xs = poss[ok]

        #pl.text(2000,i/2.+1,"b%2.2i: %1.3f" % (i,np.std(m)), size=8)


pl.arrow(0,-1,0.5*scale,0,lw=2)
pl.text(0,-2,'Half pixel')

pl.ylim([48,-3])
pl.xlim([-50,2500])

for pos in poss:
        pl.axvline(pos,color='black',ls='-.',lw=.5)

pl.xlabel("Pixel (Wavelength Direction)")
pl.ylabel("Slit #")
pl.title("Deviation of bar from predicted position")


ms = []
for i in range(nbars):
        v = means[i,:]
        v = v[np.isfinite(v)]

        ms.extend(v - np.median(v))



s = "Position"

for b in range(1, nbars, 2):
        s += "   Bar %2.0i       Bar %2.0i       " % (b, b+1)
 
s += "\n"

for p in range(npos):
        s += "%8.0f  " % (poss[p])
        for b in range(1,nbars,2):
                if np.isfinite(means[b][p]) and np.isfinite(means[b-1][p-1]) and np.isfinite(stds[b-1][p-1]) and np.isfinite(stds[b][p]):
                        s += "% 1.2f+-%1.2f  % 1.2f+-%1.2f     " % (means[b-1][p-1], stds[b-1][p-1], means[b][p], stds[b][p])
                else:
                        s += "                             "


        s += "\n"

f = open("/users/npk/desktop/for_ccs.txt","w")
f.write(s)
f.close()


