'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class expweibullpoisson(Distribution): #exponentiated weibull-poisson
    @staticmethod
    def random(a,b,g,t):
        u=ds.rg0()
        return 1/b*math.pow(-math.log(1-math.pow(1/t*math.log(u*(math.exp(t)-1)+1),1/a)),1/g)
    @staticmethod
    def pdf(a,b,g,t,x):
        return a*g*t*math.pow(b,g)*math.pow(x,g-1)/(math.exp(t)-1)*math.exp(-math.pow(b*x,g))*math.pow(1-math.exp(-math.pow(b*x,g)),a-1)*math.exp(t*math.pow(1-math.exp(-math.pow(b*x,g)),a))
    @staticmethod
    def cdf(a,b,g,t,x):
        return (math.exp(t*math.pow(1-math.exp(-math.pow(b*x,g)),a))-1)/(math.exp(t)-1)
    @staticmethod
    def median(a,b,g,t):
        return 1/b*math.pow(-math.log(1-math.pow(1/t*math.log((math.exp(t)-1)/2+1),1/a)),1/g)
    @staticmethod
    def ppf(a,b,g,t,q):
        return 1/b*math.pow(-math.log(1-math.pow(1/t*math.log(q*(math.exp(t)-1)+1),1/a)),1/g)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expweibullpoisson.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'g':ret[2],'t':ret[3]}