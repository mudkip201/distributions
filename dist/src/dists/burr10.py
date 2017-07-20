'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class burr10(Distribution):
    @staticmethod
    def random(r,k,l,s):
        return math.log(math.pow((2/(1-ds.rg0())-2)/k+1,1/r)-1)*s+l
    @staticmethod
    def pdf(r,k,l,s,x):
        xx=(x-l)/s
        return 1/s*2*math.exp(xx)*math.pow(1+math.exp(xx),r-1)*k*r/math.pow(2+(-1+math.pow(1+math.exp(xx),r))*k,2)
    @staticmethod
    def cdf(r,k,l,s,x):
        return 1-2/(2+k*math.pow(1+math.exp((x-l)/s),r)-1)
    @staticmethod
    def median(r,k,l,s):
        return math.log(math.pow(2/k+1,1/r)-1)*s+l
    @staticmethod
    def ppf(r,k,l,s,q):
        return math.log(math.pow((2/(1-q)-2)/k+1,1/r)-1)*s+l
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=burr10.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,0,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,10),(0.01,10),(-25,25),(0.01,10)]).x.tolist()
        return {'r':ret[0],'k':ret[1],'l':ret[2],'s':ret[3]}