'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class transmutedexponentiatedfrechet(Distribution):
    @staticmethod
    def random(a,b,l,t):
        return t*math.pow(-math.log(1-math.pow(((l-1)+math.sqrt((l+1)**2)-4*l*ds.rg0())/2*l,1/a)),-1/b)
    @staticmethod
    def pdf(a,b,l,t,x):
        qx=1-math.exp(-math.pow(t/x,b))
        return a*b*math.pow(t,b)*math.pow(x,-1-b)*math.exp(-math.pow(t/x,b))*math.pow(qx,a-1)*((1-l)+2*l*math.pow(qx,a))
    @staticmethod
    def cdf(a,b,l,t,x):
        return (1-math.pow(1-math.exp(-math.pow(t/x,b)),a))*(1+l*math.pow(1-math.exp(-math.pow(t/x,b)),a))
    @staticmethod
    def median(a,b,l,t):
        return t*math.pow(-math.log(1-math.pow(((l-1)+math.sqrt(l+1**2))/2*l,1/a)),-1/b)
    @staticmethod
    def ppf(a,b,l,t,q):
        return t*math.pow(-math.log(1-math.pow(((l-1)+math.sqrt((l+1)**2)-4*l*q)/2*l,1/a)),-1/b)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=transmutedexponentiatedfrechet.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2],'t':ret[3]}