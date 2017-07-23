'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class wakeby(Distribution): #if you /really/ want the kurtosis, use Mathematica/WolframAlpha
    @staticmethod
    def random(a,b,c,d,m):
        n=r.random()
        #return -a*pow(1-n,b)+c*pow(1-n,-d)+e
        return m+a*(1-math.pow(1-n,b))/b-c*(1-math.pow(1-n,-d))/d
    @staticmethod
    def mean(a,b,c,d,m):
        if(d<1):
            return a/(b+1)+c/(1-d)+m
    @staticmethod
    def median(a,b,c,d,m):
        return m+a*(1-math.pow(1/2,b))/b-c*(1-math.pow(1/2,-d))/d
    @staticmethod
    def variance(a,b,c,d,m):
        if(d<1/2):
            return a**2/((b+1)**2*(2*b+1))-2*a*c/((b+1)*(d-1)*(b-d+1)-c**2/((d-1)**2*(2*d-1)))
        return None
    @staticmethod
    def stddev(a,b,c,d,m):
        if(d<1/2):
            return wakeby.variance(a,b,c,d,m)
        return None
    @staticmethod
    def skewness(a,b,c,d,m):
        if(d<1/3):
            b1=b+1
            d1=d-1
            bd1=b-d+1
            bd2=b-2*d+1
            dd1=2*d-1
            bb1=2*b+1
            return (-2*a**3*(b-1)/(b1**3*(b*(6*b+5)+1))+(6*a**2*c*(b**2-b*(d+1)-1))/(b1**2*bb1*d1*bd1*(2*b-d+1))+(6*a*c**2*(-b*d+d**2+d-1)/(b1*d1**2*dd1*bd2*bd1))-(2*c**3*(d+1))/(d1**3*(d*(6*d-5)+1)))/math.pow(a**2/(b1**2*bb1)+2*a*c/(-b1*d1*bd1)-c**2/(d1**2*dd1),3/2)
        return None
    @staticmethod
    def ppf(a,b,c,d,m,q):
        return m+a*(1-math.pow(1-q,b))/b-c*(1-math.pow(1-q,-d))/d