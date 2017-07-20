'''
Created on Jul 16, 2017

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

class negbin(Distribution):
    @staticmethod
    def random(rr,p):
        if(rr%1!=0 or rr<=0):
            raise ValueError("rr must be a positive integer")
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        num_failure=0
        num_success=0
        while(num_failure<rr):
            if(r.random()<=p):
                num_success+=1
            else:
                num_failure+=1
        return num_success
    @staticmethod
    def kurtosis(r,p):
        if(r%1!=0 or r<=0):
            raise ValueError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        return 6/r+(1-p)**2/(p*r)
    @staticmethod
    def mean(r,p):
        if(r%1!=0 or r<=0):
            raise ValueError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        return p*r/(1-p)
    @staticmethod
    def mode(r,p):
        if(r%1!=0 or r<=0):
            raise ValueError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        if(r>1):
            return math.floor(p*(r-1)/(1-p))
        return 0
    @staticmethod
    def variance(r,p):
        if(r%1!=0 or r<=0):
            raise ValueError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        return p*r/((1-p)**2)
    @staticmethod
    def stddev(r,p):
        if(r%1!=0 or r<=0):
            raise ValueError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        return math.sqrt(p*r)/(1-p)
    @staticmethod
    def skewness(r,p):
        if(r%1!=0 or r<=0):
            raise ValueError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise ValueError("p must be between 0 and 1 exclusive")
        return (1+p)/math.sqrt(p*r)