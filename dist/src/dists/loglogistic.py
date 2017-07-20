'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import scipy.optimize as op

class loglogistic(Distribution):
    @staticmethod
    def random(aa,bb):
        n=r.random()
        while(n==1):
            n=r.random()
        return aa*math.pow((n/(1-n)),1/bb)
    @staticmethod
    def pdf(aa,bb,x):
        return (bb/aa)*math.pow(x/aa,bb-1)/math.pow(1+math.pow(x/aa,bb),2)
    @staticmethod
    def cdf(aa,bb,x):
        return 1/(1+math.pow(x/aa,-bb))
    @staticmethod
    def mean(aa,bb):
        if(bb>1):
            return aa*math.pi/bb/(math.sin(math.pi/bb))
        return None
    @staticmethod
    def median(aa,bb):
        return aa
    @staticmethod
    def mode(aa,bb):
        if(bb>1):
            return aa*math.pow((bb-1)/(bb+1),1/bb)
        return 0
    @staticmethod
    def variance(aa,bb):
        if(bb>2):
            b=math.pi/bb
            return aa**2*(2*b/math.sin(2*b)-b**2/(math.sin(b)**2))
    @staticmethod
    def stddev(aa,bb):
        if(bb>2):
            b=math.pi/bb
            return aa*math.sqrt(2*b/math.sin(2*b)-b**2/(math.sin(b)**2))
    @staticmethod
    def ppf(aa,bb,q):
        return aa*math.pow((q/(1-q)),1/bb)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=loglogistic.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0],'bb':ret[1]}