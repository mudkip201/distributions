'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp
import scipy.optimize as op

class benktanderweibull(Distribution):
    @staticmethod
    def pdf(a,b,x):
        if(b<=0 or b>1 or x<1):
            raise ValueError("b must be greater than 0 and less than or equal to 1, and x must be greater than or equal to 1")
        return math.exp(a*(1-math.pow(x,b))/b)*math.pow(x,b-1)*(1-b+a*math.pow(x,b))
    @staticmethod
    def cdf(a,b,x):
        if(b<=0 or b>1 or x<1):
            raise ValueError("b must be greater than 0 and less than or equal to 1, and x must be greater than or equal to 1")
        return 1-math.exp(a*(1-math.pow(x,b))/b)*math.pow(x,b-1)
    #@staticmethod
        #def random(a,b): #something's wrong here
        #if(b<=0 or b>1):
        #   raise ValueError("b must be greater than 0 and less than or equal to 1")
        #n=r.random()
        #return math.pow((b-1)*sp.lambertw(-(a*math.exp(-a/(b+1))*math.pow(1-n,b/(b+1)))/(b-1))/a,1/b)
    @staticmethod
    def mean(a,b):
        if(b<=0 or b>1):
            raise ValueError("b must be greater than 0 and less than or equal to 1")
        return 1+1/a
    @staticmethod
    def ppf(a,b,q):
        if(b<=0 or b>1):
            raise ValueError("b must be greater than 0 and less than or equal to 1")
        return math.pow((b-1)*sp.lambertw(-(a*math.exp(-a/(b+1))*math.pow(1-q,b/(b+1)))/(b-1))/a,1/b)
    @staticmethod
    def mle(x): #doesn't work
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=benktanderweibull.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(0.5,0.3),method='TNC',bounds=[(0.1,10),(0.1,0.999)]).x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.1,10),(0.1,0.999)]).x.tolist()
        return {'a':ret[0],'b':ret[1]}