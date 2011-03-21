
import numpy as np
import pylab as pl
import scipy as sp
import scipy.io
import os

path = "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" 
proto = "m11031%1.1i_%4.4i.fits.sav.mat"

rs = []
#for i in range(232,392):
#for i in range(232, 390):
for i in range(395,554):
        if i == 426: continue
        if os.path.exists(path + proto % (2, i)):
                ld = scipy.io.loadmat(path + proto % (2, i))
        elif os.path.exists(path + proto % (3, i)):
                ld = scipy.io.loadmat(path + proto % (3, i))
        else:
                print i
                raise Exception("No luck")

        rs.append(ld)


ixs = range(len(rs))

x = []
for r in rs:
        x.append(np.round(r["poss"]/10.)*10)

x = np.array(x)

xrs = np.unique(x[np.isfinite(x)])

bars = rs[0]["bars"]

pl.ion()
pl.figure(1)
pl.clf()
plots = []
names=[]
for i in range(len(bars)):

        if i % 2 == 0: continue
        names.append("b%2.0i" % bars[i])
        x = []
        y = []
        for r in rs:
                if r["deltas"][i] > -1: continue
                x.append(r["poss"][i])
                y.append(r["deltas"][i]+i)

        x = np.array(x).ravel()
        y = np.array(y).ravel()
        pl.plot(x,y,'-o')


#for i in range(len(bars)):
for i in [0]:
        if i % 2 == 1: continue

        yrs = xrs.copy()
        yrss = xrs.copy()
        try: finalys
        except: finalys = [[] for i in range(len(xrs))]

        for i in range(len(xrs)):
                xr = xrs[i]
                r = np.abs(x-xr) < 12
                if r.any() == False: continue
                yrs[i] = y[r].mean()
                yrss[i] = y[r].std()
                finalys[i].append(y[r].mean())
        
        pl.errorbar(xrs, yrs, yrss, fmt='^-')

        print "%3.2f" % np.std(y[0:-1])

        pf = np.poly1d(np.polyfit(x,y,1))
        print "around line: %3.2f" % np.std(y - pf(x))
        print "of mean: %3.2f" % np.std(yrs)
        try:
                print "final fit: %3.2f" % np.std(yrs[3:-3]-final_ff(xrs[3:-3]))
        except: pass



if False:
        fys = []
        sfys = []
        for fy in finalys:
                fys.append(np.mean(fy))
                sfys.append(np.std(fy))

        #pl.errorbar(xrs, fys, sfys, lw=4, fmt='k-')
        final_ff = np.poly1d(np.polyfit(xrs[3:-3],fys[3:-3],1))
        pl.plot(xrs,final_ff(xrs),'k-',lw=4)



        yl = pl.ylim()
        pl.fill([0,500,500,0,0],[yl[0],yl[0],yl[1],yl[1],yl[0]],alpha=.1,color='r')
        pl.fill([1500,2000,2000,1500,1500],[yl[0],yl[0],yl[1],yl[1],yl[0]],alpha=.1,color='r')

obars = range(1,13,2)
ebars = range(2,13,2)

vs = ["deltas", "poss"]
for p in ["o", "e"]:
        for var in vs:
            exec("%s%s = []" % (p, var))
for i in range(len(obars)):
        for r in rs:
                assert(bars[ebars[i]] == r["bars"][ebars[i]])
                assert(bars[obars[i]] == r["bars"][obars[i]])
                assert(ebars[i]-1 == obars[i])
                for p in ["o", "e"]:
                        for v in vs:
                                exec('%s%s.append(r["%s"][%sbars])' % (p, v, v, p))

print "here"
[odeltas, oposs, edeltas, eposs] = map(np.array, [odeltas, oposs, edeltas, eposs])
[odeltas, oposs, edeltas, eposs] = map(np.ravel, [odeltas, oposs, edeltas, eposs])

oddbars = obars * int(len(odeltas)/len(obars))
evenbars = ebars * int(len(odeltas)/len(obars))
