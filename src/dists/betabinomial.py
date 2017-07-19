'''
Created on Jul 15, 2017

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
import dists.beta.beta as beta
import dists.binomial.binomial as binomial

class betabinomial(Distribution):
    @staticmethod
    def random(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise ValueError("aa and bb must be bigger than 0 and n must be a positive integer")
        return binomial.random(n,beta.random(aa,bb))
    @staticmethod
    def kurtosis(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise ValueError("aa and bb must be bigger than 0 and n must be a positive integer")
        u=(aa+bb)**2*(1+aa+bb)/(n*aa*bb*(aa+bb+2)*(aa+bb+3)*(aa+bb+n))
        v=(aa+bb)*(aa+bb-1+6*n)+3*aa*bb*(n-2)+6*n**2-3*aa*bb*n*(6-n)/(aa+bb)-18*aa*bb*n**2/((aa+bb)**2)
        return u*v
    @staticmethod
    def mean(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise ValueError("aa and bb must be bigger than 0 and n must be a positive integer")
        return n*aa/(aa+bb)
    @staticmethod
    def variance(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise ValueError("aa and bb must be bigger than 0 and n must be a positive integer")
        return n*aa*bb*(aa+bb+n)/((aa+bb)**2*(aa+bb+1))
    @staticmethod
    def stddev(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise ValueError("aa and bb must be bigger than 0 and n must be a positive integer")
        return math.sqrt(n*aa*bb*(aa+bb+n)/((aa+bb)**2*(aa+bb+1)))
    @staticmethod
    def skewness(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise ValueError("aa and bb must be bigger than 0 and n must be a positive integer")
        return (aa+bb+2*n)*(bb-aa)/(aa+bb+2)*math.sqrt((1+aa+bb)/(n*aa*bb)*(n+aa+bb))