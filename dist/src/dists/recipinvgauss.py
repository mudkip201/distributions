'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.stats as st
import scipy.optimize as op

class recipinvgauss(Distribution): #reciprocal-inverse gaussian
    @staticmethod
    def random(mu):
        return st.recipinvgauss.rvs(mu)
    @staticmethod
    def pdf(mu,x):
        return 1/math.sqrt(2*math.pi*x)*math.exp(-(1-mu*x)**2/(2*x*mu**2))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=recipinvgauss.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0]}