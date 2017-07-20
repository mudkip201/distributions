'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class logisticweibull(Distribution):
    @staticmethod
    def random(a,b,l):
        return math.pow((math.pow(1/ds.rg0(),-1/l)-1)/b,1/a)
    @staticmethod
    def pdf(a,b,l,x):
        return l*a/math.pow(b,l)*math.pow(x,-l*a-1)*1/math.pow(1+math.pow(b*math.pow(x,a),-l),2)
    @staticmethod
    def cdf(a,b,l,x):
        return 1/(1+math.pow(b*math.pow(x,a),-l))
    @staticmethod
    def median(a,b,l):
        return math.pow((math.pow(2,-1/l)-1)/b,1/a)
    @staticmethod
    def ppf(a,b,l,q):
        return math.pow((math.pow(1/q,-1/l)-1)/b,1/a)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logisticweibull.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2]}