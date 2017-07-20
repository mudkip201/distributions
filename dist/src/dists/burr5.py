'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class burr5(Distribution):
    @staticmethod
    def random(r,k,l,s):
        return s*math.atan(-math.log((math.pow(ds.rg0(),-1/r)-1)/k))+l
    @staticmethod
    def pdf(r,k,l,s,x): #something's wrong here
        return r*k*math.pow(1+k/math.exp(math.tan((x-l)/s)),-r-1)*(1/math.cos((x-l)/s)**2)/math.exp(math.tan(x))/s
    @staticmethod
    def cdf(r,k,l,s,x):
        return math.pow(1+k*math.exp(-math.tan((x-l)/s)),-r)
    @staticmethod
    def median(r,k,l,s):
        return s*math.atan(-math.log((math.pow(1/2,-1/r)-1)/k))+l
    @staticmethod
    def ppf(r,k,l,s,q):
        return s*math.atan(-math.log((math.pow(q,-1/r)-1)/k))+l
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=burr5.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,0,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,10),(0.01,10),(-25,25),(0.01,math.pi)]).x.tolist()
        return {'r':ret[0],'k':ret[1],'l':ret[2],'s':ret[3]}