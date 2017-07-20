'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class expgenextremevalue(Distribution): #exponential generalized extreme value
    @staticmethod
    def random(b,k):
        if(k==0):
            return -math.log(-math.log(ds.rg0())/b)
        return (1-math.pow(-math.log(ds.rg0())/b,k))/k
    @staticmethod
    def pdf(b,k,x):
        if(k==0):
            return b*math.exp(-b*math.exp(-x))*math.exp(-x)
        return b*math.pow(1-k*x,1/k-1)*math.exp(-b*math.sqrt(1-k*x))
    @staticmethod
    def cdf(b,k,x):
        if(k==0):
            return math.exp(-b*math.exp(-x))
        return math.exp(-b*math.pow(1-k*x,1/k))
    @staticmethod
    def mean(b,k):
        if(k!=0):
            return 1/k*(1-math.pow(b,-k)*math.gamma(k+1))
    @staticmethod
    def median(b,k):
        if(k==0):
            return -math.log(-math.log(1/2)/b)
        return (1-math.pow(-math.log(1/2)/b,k))/k
    @staticmethod
    def ppf(b,k,q):
        if(k==0):
            return -math.log(-math.log(q)/b)
        return (1-math.pow(-math.log(q)/b,k))/k
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expgenextremevalue.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.1,10),(0.01,(1/max(x))-0.01)]).x.tolist()
        return {'b':ret[0],'k':ret[1]}