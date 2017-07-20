'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class arctan(Distribution):
    @staticmethod
    def pdf(lmbda,phi,x):
        if(lmbda<=0):
            raise ValueError("lambda must be greater than 0")
        return lmbda/((math.atan(lmbda*phi)+1/(2*math.pi))*(1+lmbda**2*((x-phi)**2)))
    @staticmethod
    def cdf(lmbda,phi,x):
        if(lmbda<=0):
            raise ValueError("lambda must be greater than 0")
        return 2*((math.atan(lmbda*phi)-math.atan(-x*lmbda+lmbda*phi))/(2*math.atan(lmbda*phi)+math.pi))
    @staticmethod
    def random(lmbda,phi):
        if(lmbda<=0):
            raise ValueError("lambda must be greater than 0")
        n=ds.rg0()
        return (lmbda*phi+math.tan(-math.atan(lmbda*phi)+n*math.atan(lmbda*phi)+n*math.pi/2))/lmbda
    @staticmethod
    def kurtosis(lmbda,phi):
        if(lmbda<=0):
            raise ValueError("lambda must be greater than 0")
        return None
    @staticmethod
    def mean(lmbda,phi):
        if(lmbda<=0):
            raise ValueError("lambda must be greater than 0")
        return None
    @staticmethod
    def variance(lmbda,phi):
        if(lmbda<=0):
            raise ValueError("lambda must be greater than 0")
        return float("infinity")
    @staticmethod
    def stddev(lmbda,phi):
        if(lmbda<=0):
            raise ValueError("lambda must be greater than 0")
        return float("infinity")
    @staticmethod
    def ppf(lmbda,phi,q):
        if(lmbda<=0):
            raise ValueError("lambda must be greater than 0")
        return (lmbda*phi+math.tan(-math.atan(lmbda*phi)+q*math.atan(lmbda*phi)+q*math.pi/2))/lmbda
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=arctan.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='L-BFGS-B',bounds=[(0,1),(-10000,10000)]).x.tolist()
        return {'lambda':ret[0],'phi':ret[1]}