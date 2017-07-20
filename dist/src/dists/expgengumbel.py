'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class expgengumbel(Distribution): #exponentiated generalized gumbel
    @staticmethod
    def random(a,b,m,s):
        return m-s*math.log(-math.log(1-math.pow(1-math.pow(ds.rg0(),1/b),1/a)))
    @staticmethod
    def pdf(a,b,m,s,x):
        return math.pow(1-math.pow(1-math.exp(-math.exp(-(x-m)/s)),a),b)
    @staticmethod
    def cdf(a,b,m,s,x):
        return a*b/s*math.exp(-((x-m)/s+math.exp(-(x-m)/s)))*math.pow(1-math.exp(-math.exp(-(x-m)/s)),a-1)*math.pow(1-math.pow(1-math.exp(-math.exp(-(x-m)/s)),a),b-1)
    @staticmethod
    def median(a,b,m,s):
        return m-s*math.log(-math.log(1-math.pow(1-math.pow(1/2,1/b),1/a)))
    @staticmethod
    def ppf(a,b,m,s,q):
        return m-s*math.log(-math.log(1-math.pow(1-math.pow(q,1/b),1/a)))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expgengumbel.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,0,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'m':ret[2],'s':ret[3]}