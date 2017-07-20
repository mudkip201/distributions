'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp
import scipy.optimize as op
import dists.f.f as f

class fisherz(Distribution):
    @staticmethod
    def random(n,m):
        return math.log(f.random(n,m))/2
    @staticmethod
    def pdf(n,m,z):
        return 2*math.pow(n,n/2)*math.pow(m,m/2)/sp.beta(n/2,m/2)*math.exp(n*z)/math.pow(n*math.exp(2*z)+m,(n+m)/2)
    @staticmethod
    def mean(n,m):
        return m/(m-2)
    @staticmethod
    def mode(n,m):
        return m/(m+2)*(n-2)/n
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=fisherz.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'n':ret[0],'m':ret[1]}