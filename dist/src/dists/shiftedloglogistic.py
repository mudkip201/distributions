'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class shiftedloglogistic(Distribution):
    @staticmethod
    def random(xi,mu,sigma):
        n=ds.rg0()
        return (math.pow(1/n-2,-xi)*(xi*mu*math.pow(1/n-2,xi)+sigma))/xi
    @staticmethod
    def pdf(xi,mu,sigma,x):
        z=(x-mu)/sigma
        return math.pow(1+xi*z,-1/xi-1)/(sigma*(1+math.pow(1+xi*z,-1/xi))**2)
    @staticmethod
    def cdf(xi,mu,sigma,x):
        z=(x-mu)/sigma
        return 1/(1+math.pow(1+xi*z,-1/xi))
    @staticmethod
    def mean(xi,mu,sigma):
        return mu+sigma/xi*(math.pi*xi/math.sin(math.pi*xi)-1)
    @staticmethod
    def median(xi,mu,sigma):
        return mu
    @staticmethod
    def mode(xi,mu,sigma):
        return mu+sigma/xi*(math.pow((1-xi)/(1+xi),xi)-1)
    @staticmethod
    def variance(xi,mu,sigma):
        aa=math.pi*xi
        return sigma**2/xi**2*(2*aa/math.sin(2*aa)-(aa/math.sin(aa))**2)
    @staticmethod
    def stddev(xi,mu,sigma):
        aa=math.pi*xi
        return sigma/xi*math.sqrt(2*aa/math.sin(2*aa)-(aa/math.sin(aa))**2)
    @staticmethod
    def ppf(xi,mu,sigma,q):
        return (math.pow(1/q-2,-xi)*(xi*mu*math.pow(1/q-2,xi)+sigma))/xi
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=shiftedloglogistic.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'xi':ret[0],'mu':ret[1],'sigma':ret[2]}
