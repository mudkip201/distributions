'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.optimize as op
import dists.gamma.gamma as gamma

class invgamma(Distribution): #inverse gamma (Pearson V)
    @staticmethod
    def random(a,b):
        return 1/gamma.random(a,b)
    @staticmethod
    def pdf(a,b,x):
        return math.pow(b,a)/math.gamma(a)*math.pow(x,-a-1)*math.exp(-b/x)
    @staticmethod
    def kurtosis(a,b):
        if(a>4):
            return (30*a-66)/((a-3)*(a-4))
    @staticmethod
    def mean(a,b):
        if(a>1):
            return b/(a-1)
    @staticmethod
    def mode(a,b):
        return b/(a+1)
    @staticmethod
    def variance(a,b):
        if(a>2):
            return b**2/((a-1)**2*(a-2))
    @staticmethod
    def stddev(a,b):
        if(a>2):
            return math.sqrt(b**2/((a-1)**2*(a-2)))
    @staticmethod
    def skewness(a,b):
        if(a>3):
            return 4*math.sqrt(a-2)/(a-3)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=invgamma.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1]}
