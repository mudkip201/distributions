'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import scipy.optimize as op

class kumaraswamy(Distribution):
    @staticmethod
    def pdf(a,b,x):
        if(a<=0 or b<=0):
            raise ValueError("a and b must be positive")
        if(x<0 or x>1):
            raise ValueError("x must be between 0 and 1 inclusive")
        return a*b*math.pow(x,a-1)*math.pow(1-math.pow(x,a),b-1)
    @staticmethod
    def cdf(a,b,x):
        if(a<=0 or b<=0):
            raise ValueError("a and b must be positive")
        if(x<0 or x>1):
            raise ValueError("x must be between 0 and 1 inclusive")
        return 1-math.pow(1-math.pow(x,a),b)
    @staticmethod
    def random(a,b):
        if(a<=0 or b<=0):
            raise ValueError("a and b must be positive")
        return math.pow(1-math.pow(1-r.random(),1/b),1/a)
    @staticmethod
    def mean(a,b):
        if(a<=0 or b<=0):
            raise ValueError("a and b must be positive")
        return b*math.gamma(1+1/a)*math.gamma(b)/math.gamma(1+1/a+b)
    @staticmethod
    def median(a,b):
        if(a<=0 or b<=0):
            raise ValueError("a and b must be positive")
        return math.pow(1-math.pow(2,-1/b),1/a)
    @staticmethod
    def mode(a,b):
        if(a<=0 or b<=0):
            raise ValueError("a and b must be positive")
        if(a>=1 and b>=1 and (a!=1 or b!=1)):
            return math.pow((a-1)/(a*b-1),1/a)
        return None
    @staticmethod
    def variance(a,b):
        if(a<=0 or b<=0):
            raise ValueError("a and b must be positive")
        return (b*math.gamma(1+2/a)*math.gamma(b)/math.gamma(1+b+2/a))-(b*math.gamma(1+1/a)*math.gamma(b)/math.gamma(1+b+1/a))**2
    @staticmethod
    def stddev(a,b):
        if(a<=0 or b<=0):
            raise ValueError("a and b must be positive")
    @staticmethod
    def skewness(a,b):
        return b*math.gamma(b)*(a**3*math.gamma((a+3)/a)/math.gamma(b+3/a+1)-6*a*math.gamma(1/a)*math.gamma(2/a)*math.gamma(b+1)/(math.gamma(b+1/a+1)*math.gamma(b+2/a+1))+2*math.gamma(1/a)**3*math.gamma(b+1)**2/math.gamma(b+1/a+1)**3)/(a**3*math.pow(math.gamma(1+2/a)*math.gamma(b+1)/math.gamma(b+2/a+1)-math.gamma(1+1/a)**2*math.gamma(b+1)**2/math.gamma(b+1/a+1)**2,3/2))
    @staticmethod
    def kurtosis(a,b):
        gb4=math.gamma(b+4/a+1)
        gb3=math.gamma(b+3/a+1)
        gb2=math.gamma(b+2/a+1)
        gb1=math.gamma(b+1/a+1)
        gb=math.gamma(b+1)
        ga1=math.gamma(1/a)
        ga2=math.gamma(2/a)
        ga3=math.gamma(3/a)
        k=(gb2*(a**4+math.gamma((a+4)/a)*gb1**4*gb2*gb3-3*ga1*gb*(4*a**2*ga3*gb2*gb1**3+ga1*gb*(ga1**2*gb*gb2-4*a*ga2*gb1**2)*gb3)*gb4))
        k/=gb*(a**4*math.gamma((a+2)/a)**2*gb1**4+ga1**2*gb*gb2*(ga1**2*gb*gb2-4*a*ga2*gb1**2))*gb3*gb4
    @staticmethod
    def ppf(a,b,q):
        if(a<=0 or b<=0):
            raise ValueError("a and b must be positive")
        return math.pow(1-math.pow(1-q,1/b),1/a)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamy.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1]}