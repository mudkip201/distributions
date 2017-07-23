'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class transmutedcompweibullgeo(Distribution): #transmuted complementary weibull geometric
    @staticmethod
    def random(a,b,d,g):
        q=ds.rg0()
        return math.pow(math.log((2*a*a-2*a*q*(a-1)-a*(1+d)+a*math.sqrt(1+d*(d-4*q+2)))/(2*a*a*(1-q))),1/b)/g
    @staticmethod
    def pdf(a,b,d,g,x):
        return a*b*g*math.pow(g*x,b-1)*math.exp(-math.pow(g*x,b))*(a*(1-d)-(a-a*d-d-1)*math.exp(-math.pow(g*x,b)))/math.pow(a+(1-a)*math.exp(-math.pow(g*x,b)),3)
    @staticmethod
    def median(a,b,d,g):
        return math.pow(math.log((2*a**2-a*(a-1)-a*(1+d)+a*math.sqrt(1+d**2))/(a**2)),1/b)/g
    @staticmethod
    def ppf(a,b,d,g,q):
        return math.pow(math.log((2*a*a-2*a*q*(a-1)-a*(1+d)+a*math.sqrt(1+d*(d-4*q+2)))/(2*a*a*(1-q))),1/b)/g
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=transmutedcompweibullgeo.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'d':ret[2],'g':ret[3]}
