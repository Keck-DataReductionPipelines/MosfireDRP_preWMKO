''' Test of beta mask-level fitting '''

import numpy as np
import pylab as pl

from MOSFIRE import Fit, Wavelength

reload(Wavelength)
reload(Fit)

sols = np.load('/Users/npk/desktop/c9_reduce/npk_calib3_q1700_pa_0/lambda_coeffs_m110323_2718.npy')


filter_fun = (lambda x:
        (x[1] is not None) and
        (x[1][0] < 2e-5) and
        (x[2] < .1) and
        (x[3] == True))


betas = []
pixs = []

all_betas = []
all_pixs = []

for solution in sols:
    sol_2d = solution["2d"]
    print "Slit {0}".format(solution["slitno"])

    ff = filter(filter_fun, sol_2d)
    ar = np.array(map(lambda x: x[0], ff))

    p = ar[:,5]
    b = ar[:,1]
    ss = np.argsort(p)
    pixs.append(p[ss])
    betas.append(b[ss])


a = pixs[0]
b = pixs[1]
del pixs[0]
del pixs[1]
pixs.insert(0, np.append(a,b))

a = betas[0]
b = betas[1]
del betas[0]
del betas[1]
betas.insert(0, np.append(a,b))

ind = np.argsort(pixs[0])
pixs[0] = pixs[0][ind]
betas[0] = betas[0][ind]

all_pixs = np.concatenate(pixs)
all_betas = np.concatenate(betas)

pars = [1500, 20, 100000, 10000]
pars.extend(np.zeros(len(pixs)))
pars = np.array(pars)

parinfo = []
for par in pars:
    parinfo.append({"fixed": 0, "value": par})

pl.ion()
y = Wavelength.beta_model(pars, pixs)

parinfo[0]["fixed"] = 1
parinfo[1]["fixed"] = 1
parinfo[2]["fixed"] = 1
parinfo[3]["fixed"] = 1
merit_fun = Fit.mpfit_residuals(Wavelength.beta_model)
lsf = Fit.mpfit_do(merit_fun, pixs, all_betas, parinfo)
parinfo = []
for param in lsf.params:
    parinfo.append({"fixed": 0, "value": param})
lsf = Fit.mpfit_do(merit_fun, pixs, all_betas, parinfo)
print lsf


pl.figure(2)
pl.clf()
pl.plot(all_pixs, all_betas, '.')
y = Wavelength.beta_model(lsf.params, pixs)
d = np.abs(np.diff(y))
rois = np.where(d>0.01)[0]
prev = 0
for roi in rois:
    pl.plot(all_pixs[prev+1:roi-1], y[prev+1:roi-1])
    prev=roi

pl.figure(3)
pl.clf()
y = Wavelength.beta_model(lsf.params, pixs)
pl.plot(all_pixs, all_betas-y)

