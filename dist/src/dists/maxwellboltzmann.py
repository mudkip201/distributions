'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class maxwellboltzmann(Distribution):
    @staticmethod
    def pdf(aa,x):
        if(aa<=0 or x<=0):
            raise ValueError("all inputs must be positive")
        return math.sqrt(2/math.pi)*(x**2*math.exp(x**2/(2*aa**2)))/(aa**3)
    @staticmethod
    def cdf(aa,x):
        if(aa<=0 or x<=0):
            raise ValueError("all inputs must be positive")
        return math.erf(x/(math.sqrt(2)*aa))-math.sqrt(2/math.pi)*(x*math.exp(x**2/(2*aa**2)))/aa
    @staticmethod
    def random(aa):
        if(aa<=0):
            raise ValueError("aa must be positive")
        r=(-math.log(ds.rg0()))**2
        w1=(ds.rg0())**2
        w2=(ds.rg0())**2
        while(w1+w2>1):
            w1=(ds.rg0())**2
            w2=(ds.rg0())**2
        r2=r-(w1/(w1+w2))*math.log(ds.rg0())
        return aa*math.sqrt(2*r2)
    @staticmethod
    def kurtosis(aa):
        if(aa<=0):
            raise ValueError("aa must be positive")
        return 4*(-96+40*math.pi-3*(math.pi)**2)/((3*math.pi-8)**2)
    @staticmethod
    def mean(aa):
        if(aa<=0):
            raise ValueError("aa must be positive")
        return 2*aa*math.sqrt(2/math.pi)
    @staticmethod
    def mode(aa):
        if(aa<=0):
            raise ValueError("aa must be positive")
        return math.sqrt(2)*aa
    @staticmethod
    def variance(aa):
        if(aa<=0):
            raise ValueError("aa must be positive")
        return aa**2*(3*math.pi-8)/math.pi
    @staticmethod
    def stddev(aa):
        if(aa<=0):
            raise ValueError("aa must be positive")
        return aa*math.sqrt((3*math.pi-8)/math.pi)
    @staticmethod
    def entropy(aa):
        if(aa<=0):
            raise ValueError("aa must be positive")
        return math.log(aa*math.sqrt(2*math.pi))+ds.euler_gamma-1/2
    @staticmethod
    def skewness(aa):
        if(aa<=0):
            raise ValueError("aa must be positive")
        return 2*math.sqrt(2)*(16-5*math.pi)/math.pow(3*math.pi-8,3/2)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=maxwellboltzmann.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0]}