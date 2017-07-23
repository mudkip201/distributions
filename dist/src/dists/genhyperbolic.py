'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class genhyperbolic(Distribution):
    @staticmethod
    def pdf(a,b,d,l,m,x):
        g=math.sqrt(a**2-b**2)
        ff=math.pow(g/d,l)*math.exp(b*(x-m))/(math.sqrt(2*math.pi)*sp.kv(l,d*g))
        ff*=sp.kv(l-1/2,a*math.sqrt(d**2+(x-m)**2))/math.pow(math.sqrt(d**2+(x-m)**2)/a,1/2-l)
        return ff
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,b,d,l,m):
        g=math.sqrt(a**2-b**2)
        return m+d*b*sp.kv(l+1,d*g)/(g*sp.kv(l,d*g))
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a,b,d,l,m):
        g=math.sqrt(a**2+b**2)
        vv=d*sp.kv(l+1,d*g)/(g*sp.kv(l,d*g))
        vv+=b**2*d**2/g**2*(sp.kv(l+2,d*g)/sp.kv(l,d*g)-(sp.kv(l+1,d*g)/sp.kv(l,d*g))**2)
        return vv
    @staticmethod
    def stddev(a,b,d,l,m):
        return math.sqrt(genhyperbolic.variance(a,b,d,l,m))
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