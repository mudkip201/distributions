'''
Created on Jul 19, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import numpy as np
from numpy import random as r
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as op

class modburr3(Distribution):
    @staticmethod
    def pdf(dist,aa,bb,gg,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        g=dist.pdf(*args,**kwargs)
        H=G/(1-G)
        return aa*bb*math.pow(H,-bb+1)*math.pow(1+g*math.pow(H,-bb),-aa/gg-1)*g/G**2
    @staticmethod
    def cdf(dist,aa,bb,g,*args,**kwargs):
        return math.pow(1+g*math.pow((dist.cdf(*args,**kwargs))/(1-dist.cdf(*args,**kwargs)),-bb),-aa/g)
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass

class kumaraswamy(Distribution):
    @staticmethod
    def pdf(dist,a,b,*args,**kwargs):
        g=dist.pdf(*args,**kwargs)
        G=dist.cdf(*args,**kwargs)
        return a*b*g*math.pow(G,a-1)*math.pow(1-math.pow(G,a),b-1)
    @staticmethod
    def cdf(dist,a,b,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        return 1-math.pow(1-math.pow(G,a),b)
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass    

class weibull(Distribution):
    @staticmethod
    def pdf(dist,l,*args,**kwargs):
        g=dist.pdf(*args,**kwargs)
        G=dist.cdf(*args,**kwargs)
        return l*g*math.pow(G,l-1)/math.pow(1-G,l+1)*math.exp(-math.pow(G/(1-G),l))
    @staticmethod
    def cdf(dist,l,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        return 1-math.exp(-math.pow(G/(1-G),l))
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass
    
class beta(Distribution):
    @staticmethod
    def pdf(dist,a,b,*args,**kwargs):
        g=dist.pdf(*args,**kwargs)
        G=dist.cdf(*args,**kwargs)
        return 1/sp.beta(a,b)*g*math.pow(G,a-1)*math.pow(1-G,b-1)
    @staticmethod
    def cdf(dist,a,b,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        return sp.betainc(a,b,G)
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass
    
class hyperboliccosine(Distribution):
    @staticmethod
    def pdf(dist,a,*args,**kwargs):
        g=dist.pdf(*args,**kwargs)
        G=dist.cdf(*args,**kwargs)
        return 2*a*math.exp(a)/(math.exp(2*a)-1)*g*math.cosh(a*G)
    @staticmethod
    def cdf(dist,a,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        return 2*math.exp(a)/(math.exp(2*a)-1)*math.sinh(a*G)
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf(dist,a,*args,**kwargs):
        aargs=list(args)[:-1]
        return dist.ppf(tuple(aargs),1/a*math.asinh(math.exp(2*a-1)/(2*math.exp(a))*args[-1]),**kwargs)
        pass
    @staticmethod
    def mle():
        pass
    
class transmutedgeometric(Distribution):
    @staticmethod
    def pdf(dist,l,t,*args,**kwargs):
        g=dist.pdf(*args,**kwargs)
        G=dist.cdf(*args,**kwargs)
        return t*g/math.pow(1+(t-1)*2,2)*(1+l-2*l*t*(1-G)/(1+(t-1)*G))
    @staticmethod
    def cdf(dist,l,t,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        return t*G/(1+(t-1)*G)*(1+l(1-G)/(1+(t-1)*G))
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass
    
    
class expgengompertz(Distribution):
    @staticmethod
    def pdf(dist,a,g,t,*args,**kwargs):
        g=dist.pdf(*args,**kwargs)
        G=dist.cdf(*args,**kwargs)
        p=a*t*g*math.exp(t/g*(1-math.pow(1-G,-g)))
        p/=math.pow(1-G,1+g)
        return p*math.pow(1-math.exp(t/g*(1-math.pow(1-G,-g))),a-1)
    @staticmethod
    def cdf(dist,a,g,t,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        return math.pow(1-math.exp(t/g*(1-math.pow(1-G,-g))),a)
    @staticmethod
    def random(dist,a,g,t,*args,**kwargs):
        return dist.ppf(args,1-math.pow(1-g/t*math.log(1-math.pow(ds.rg0(),1/a))),**kwargs)
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median(dist,a,g,t,*args,**kwargs):
        return dist.ppf(args,1-math.pow(1-g/t*math.log(1-math.pow(1/2,1/a))),**kwargs)
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf(dist,a,g,t,*args,**kwargs):
        aargs=list(args)[:-1]
        return dist.ppf(tuple(aargs),1-math.pow(1-g/t*math.log(1-math.pow(args[-1],1/a))),**kwargs)
    @staticmethod
    def mle():
        pass
    
class genoddloglogistc(Distribution):
    @staticmethod
    def pdf(dist,a,t,*args,**kwargs):
        g=dist.pdf(*args,**kwargs)
        G=dist.cdf(*args,**kwargs)
        return a*t*g*math.pow(G,a*t-1)*math.pow(1-math.pow(G,t),a-1)/math.pow(math.pow(G,a*t)+math.pow(1-math.pow(G,t),a),2)
    @staticmethod
    def cdf(dist,a,t,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        return math.pow(G,a*t)/(math.pow(G,a*t)+math.pow(1-math.pow(G,t),a))
    @staticmethod
    def random(dist,a,t,*args,**kwargs):
        q=ds.rg0()
        return dist.ppf(args,q=math.pow(math.pow(q/(1-q),1/a)/(1+math.pow(q/(1-q),1/a)),1/t),**kwargs)
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median(dist,a,t,*args,**kwargs):
        return dist.ppf(args,q=math.pow(1/2,1/t),**kwargs)
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf(dist,a,t,*args,**kwargs):
        aargs=list(args)[:-1]
        q=args[0]
        return dist.ppf(tuple(aargs),q=math.pow(math.pow(q/(1-q),1/a)/(1+math.pow(q/(1-q),1/a)),1/t),**kwargs)
    @staticmethod
    def mle():
        pass
