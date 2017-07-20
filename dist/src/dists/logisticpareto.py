'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class logisticpareto(Distribution):
    @staticmethod
    def random(k,l,t):
        return t*math.pow(math.exp(math.pow(1/ds.rg0()-1,-1/l)),1/k)
    @staticmethod
    def pdf(k,l,t,x):
        return l*k/x*math.pow(math.log(math.pow(x/t,k)),-l-1)/math.pow(1+math.pow(math.log(math.pow(x/t,k)),-l),2)
    @staticmethod
    def cdf(k,l,t,x):
        return 1/(1+math.pow(math.log(math.pow(x/t,k)),-l))
    @staticmethod
    def median(k,l,t):
        return t*math.pow(math.exp(math.pow(1,-1/l)),1/k)
    @staticmethod
    def ppf(k,l,t,q):
        return t*math.pow(math.exp(math.pow(1/q-1,-1/l)),1/k)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logisticpareto.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'k':ret[0],'l':ret[1],'t':ret[2]}