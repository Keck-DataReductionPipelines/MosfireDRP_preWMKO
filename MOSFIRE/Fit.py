'''
Fitting code used by a variety of MOSFIRE applications.

Written in March 2011 by npk
'''
from scipy.special import erf
import scipy.optimize as optimize
import numpy as np
import pylab as pl

import unittest

def find_max(v):
        return np.where(v == v.max())[0][0]

def find_min(v):
        return np.where(v == v.min())[0]

def slit_edge_fun(x, s):
        ''' The edge of a slit, convolved with a Gaussian, is well fit by the error function. slit_edge_fun is a reexpression of the error function in the classica gaussian "sigma" units'''
        sq2 = np.sqrt(2)
        sig = sq2 * s

        return np.sqrt(np.pi/2.) * s * erf(x/sig) 

def fit_bar_edge(p, x):
        ''' 
        Fitting function for a bar edge
        '''

        return p[0] + np.radians(4.2) * x

def fit_single(p, x):
        '''
        The fitting function used by do_fit_single. This is a single slit edge
        
        p[0] ---> Sigma
        p[1] ---> Horizontal offset
        p[2] ---> Multipicative offset
        p[3] ---> Additive offset
        '''
        return slit_edge_fun(x - p[1], p[0]) * p[2] + p[3]

def fit_pair(p, x):
        '''
        The fitting function ussed by "do_fit". The sum of two edge functions.

        p[0] ---> Sigma
        p[1] ---> Horizontal offset
        p[2] ---> Multipicative offset
        p[3] ---> Additive offset
        p[4] ---> Width of slit
        '''
        return slit_edge_fun(x - p[1], p[0]) * p[2] + p[3] - slit_edge_fun(x - p[1] - p[4], p[0]) * p[2]

def residual(p, x, y, f):
        '''The square of residual is minimized by the least squares fit. Formally this is (f(x | p) - y)**2'''
        return f(p, x) - y

def residual_single(p, x, y):
        '''Convenience funciton around residual'''
        return residual(p, x, y, fit_single)

def residual_pair(p, x, y):
        '''Convenience funciton around residual'''
        return residual(p, x, y, fit_pair)

def residual_bar_edge(p, x, y):
        return residual(p, x, y, fit_bar_edge)

def do_fit(data, residual_fun=residual_single):
        '''do_fit estimates parameters of fit_pair or fit_single.
        
        Use as follows:

        p0 = [0.5, 6, 1.1, 3, 1]
        ys = fit_single(p0, xs)
        lsf = do_fit(ys, residual_single)
        res = np.sum((lsf[0] - p0)**2)

        '''


        xs = np.arange(len(data))

        if residual_fun==residual_single:
                if data[0] > data[-1]:
                        p0 = [0.5, len(data)/2., max(data), 0.0, 3.0]
                else:
                        p0 = [0.5, len(data)/2., -max(data), 0.0, 3.0]
        elif residual_fun==residual_pair:
                p0 = [0.5, find_max(data), max(data), 0.0, 4.0]


        lsf = optimize.leastsq(residual_fun, p0, args=(xs, data), full_output=True)

        return lsf

def do_fit_edge(xs, ys):

        p0 = [ys.mean()]

        return optimize.leastsq(residual_bar_edge, p0, args=(xs, ys))

class TestFitFunctions(unittest.TestCase):

        def setUp(self):
                pass

        def test_do_fit(self):
                import random
                sa = self.assertTrue

                xs = np.arange(19)

                p0 = [0.5, 6, 1.1, 3, 1]
                ys = fit_single(p0, xs)
                lsf = do_fit(ys, residual_single)
                res = np.sum((fit_single(lsf[0],xs) - ys)**2)
                sa(res < 0.001)

        def test_do_fit2(self):
                sa = self.assertTrue
                p0 = [0.5, 6, 1.1, 3, 3]
                xs = np.arange(15)
                ys = fit_pair(p0, xs)

                lsf = do_fit(ys, residual_pair)
                res = np.sum((lsf[0] - p0)**2)
                sa(res < 0.001)


        def test_lhs_v_rhs(self):
                sa = self.assertTrue

                p0 = [0.5, 5, 1.1,.7,0]
                pn0 = [0.5, 5, -1.1,.7,0]
                xs = np.arange(25)
                ys = fit_single(p0, xs)
                lsf = do_fit(ys, residual_single)
                print lsf[0]

                ys = fit_single(pn0, xs)
                lsf = do_fit(ys, residual_single)
                print lsf[0]

if __name__ == '__main__':
        unittest.main()

