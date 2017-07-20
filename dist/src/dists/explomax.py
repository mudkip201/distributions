'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class explomax(Distribution): #exponential lomax
    @staticmethod
    def random(a,b,l):
        return b*(math.pow(-math.log(1-ds.rg0())/l,1/a)-1)
    @staticmethod
    def pdf(a,b,l,x):
        return a*l/b*math.pow(b/(x+b),-a+1)*math.exp(-l*math.pow(b/(x+b),-a))
    @staticmethod
    def cdf(a,b,l,x):
        return 1-math.exp(-l*math.pow(b/(x+b),-a))
    @staticmethod
    def median(a,b,l):
        return b*(math.pow(-math.log(1/2)/l,1/a)-1)
    @staticmethod
    def ppf(a,b,l,q):
        return b*(math.pow(-math.log(1-q)/l,1/a)-1)
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=explomax.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2]}