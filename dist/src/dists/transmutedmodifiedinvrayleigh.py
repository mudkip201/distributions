'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class transmutedmodifiedinvrayleigh(Distribution): #transmuted modified inverse rayleigh
    @staticmethod
    def random(a,b,l):
        return 2*b/(-a+math.sqrt(a**2-4*b*math.log(((1+l)-math.sqrt((1+l)**2)-4*l*ds.rg0())/(2*l))))
    @staticmethod
    def pdf(a,b,l,x):
        return (a+2*b/x)*(1/x**2)*math.exp(-a/x-b/x**2)*(1+l-2*l*math.exp(-a/x-b/x**2))
    @staticmethod
    def cdf(a,b,l,x):
        return math.exp(-a/x-b/x**2)*(1+l-l*math.exp(-a/x-b/x**2))
    @staticmethod
    def median(a,b,l):
        return 2*b/(-a+math.sqrt(a**2-4*b*math.log(((1+l)-math.sqrt(1+l**2))/(2*l))))
    @staticmethod
    def ppf(a,b,l,q):
        return 2*b/(-a+math.sqrt(a**2-4*b*math.log(((1+l)-math.sqrt((1+l)**2)-4*l*q)/(2*l))))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=transmutedmodifiedinvrayleigh.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2]}