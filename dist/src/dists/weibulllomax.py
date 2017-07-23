'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class weibulllomax(Distribution):
    @staticmethod
    def random(a,b,aa,bb):
        return bb*(math.pow(math.pow(-math.log(1-ds.rg0())/a,1/b)+1,1/aa)-1)
    @staticmethod
    def pdf(a,b,aa,bb,x):
            return a*b*aa/bb*math.pow(1+x/bb,b*aa-1)*math.pow(1-math.pow(1+x/bb,-aa),b-1)*math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))
    @staticmethod
    def cdf(a,b,aa,bb,x):
        return 1-math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))
    @staticmethod
    def median(a,b,aa,bb):
        return bb*(math.pow(math.pow(-math.log(1/2)/a,1/b)+1,1/aa)-1)
    @staticmethod
    def ppf(a,b,aa,bb,q):
        return bb*(math.pow(math.pow(-math.log(1-q)/a,1/b)+1,1/aa)-1)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=weibulllomax.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'aa':ret[2],'bb':ret[3]}