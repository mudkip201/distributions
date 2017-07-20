'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.special as sp

class logdagum(Distribution):
    @staticmethod
    def pdf(b,d,l,x):
        return b*l*d*math.exp(-d*x)*math.pow(1+l*math.exp(-d*x),-b-1)
    @staticmethod
    def cdf(b,d,l,x):
        return math.pow(1+l*math.exp(-d*x),-b)
    @staticmethod
    def random(b,d,l):
        u=ds.rg0()
        return math.log(l/(math.pow(u,-1/b)-1))/d
    @staticmethod
    def mean(b,d,l):
        return (math.log(l)+sp.digamma(b)-sp.digamma(1))/d
    @staticmethod
    def median(b,d,l):
        return math.log(l/(math.pow(1/2,-1/b)-1))/d
    @staticmethod
    def mode(b,d,l):
        if(b*d>1):
            return math.log(l*b)/d
        return None
    @staticmethod
    def variance(b,d,l):
        return ((sp.polygamma(3,b)+sp.polygamma(3,1))+math.pow(logdagum.mean(b,d,l),2))/math.pow(d,2)-math.pow(logdagum.mean(b,d,l),2)
    @staticmethod
    def stddev(b,d,l):
        return math.sqrt(logdagum.variance(b,d,l))
    @staticmethod
    def kurtosis(b,d,l):
        e1=sp.polygamma(5,b)+sp.polygamma(5,1)
        e2=3*(sp.polygamma(3,b)+sp.polygamma(3,1))**2
        e3=4*(math.log(l)+sp.digamma(b)-sp.digamma(1))*(sp.polygamma(4,b)-sp.polygamma(4,1))
        e4=6*(math.log(l)+sp.digamma(b)-sp.digamma(1))**2*(sp.polygamma(3,b)+sp.polygamma(3,1))
        e5=(math.log(l)+sp.digamma(b)-sp.digamma(1))**4
        return (e1+e2+e3+e4+e5)/(d**4)-math.pow(logdagum.mean(b,d,l),4)
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(b,d,l):
        e1=sp.polygamma(4,b)-sp.polygamma(4,1)
        e2=math.pow(math.log(l)+sp.digamma(b)-sp.digamma(1),3)
        e3=3*(math.log(l)+sp.digamma(b)-sp.digamma(1))*(sp.polygamma(3,b)+sp.polygamma(3,1))
        return (e1+e2+e3)/(d**3)-logdagum.mean(b,d,l)**3
    @staticmethod
    def ppf(b,d,l,q):
        return math.log(l/(math.pow(q,-1/b)-1))/d
    @staticmethod
    def mle():
        pass