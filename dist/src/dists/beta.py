'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp
import dists.gamma.gamma as gamma

class beta(Distribution):
    @staticmethod
    def pdf(aa,bb,x):
        return math.pow(x,aa-1)*math.pow(1-x,bb-1)/sp.beta(aa,bb)
    @staticmethod
    def cdf(aa,bb,x):
        return sp.betainc(aa,bb,x)*sp.beta(aa,bb)
    @staticmethod
    def random(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be bigger than 0")
        x=gamma.random(aa,1)
        return x/(x+gamma.random(bb,1))
    @staticmethod
    def kurtosis(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be bigger than 0")
        return 6*((aa-bb)**2*(aa+bb+1)-aa*bb*(aa+bb+2))/(aa*bb*(aa+bb+2)*(aa+bb+3))
    @staticmethod
    def mean(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be bigger than 0")
        return aa/(aa+bb)
    @staticmethod
    def median(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be bigger than 0")
        return sp.betaincinv(aa,bb,1/2)
    @staticmethod
    def variance(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be bigger than 0")
        return aa*bb/((aa+bb)**2*(aa+bb+1))
    @staticmethod
    def stddev(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be bigger than 0")
        return math.sqrt(aa*bb/((aa+bb)**2*(aa+bb+1)))
    @staticmethod
    def skewness(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be bigger than 0")
        return 2*(bb-aa)*math.sqrt(aa+bb+1)/((aa+bb+2)*math.sqrt(aa*bb))