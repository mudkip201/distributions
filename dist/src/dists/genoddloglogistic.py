'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class genoddloglogistic(Distribution): #generalized odd log-logistic
    @staticmethod
    def random(a,t):
        u=ds.rg0()
        return math.pow(math.pow(u/(1-u),1/a)/(1+math.pow(u/(1-u),1/a)),1/t)
    @staticmethod
    def median(a,t):
        return math.pow(1/2,1/t)
    @staticmethod
    def ppf(a,t,q):
        return math.pow(math.pow(q/(1-q),1/a)/(1+math.pow(q/(1-q),1/a)),1/t)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=genoddloglogistic.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'t':ret[1]}
