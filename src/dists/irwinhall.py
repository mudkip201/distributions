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

class irwinhall(Distribution):
    @staticmethod
    def random(n):
        if(n<=0 or n%1!=0):
            raise ValueError("n must be a positive integer")
        sum_=0
        for _ in range(n):
            sum_+=r.random()
        return sum_
    @staticmethod
    def kurtosis(n):
        if(n<=0 or n%1!=0):
            raise ValueError("n must be a positive integer")
        return -6/(5*n)
    @staticmethod
    def mean(n):
        if(n<=0 or n%1!=0):
            raise ValueError("n must be a positive integer")
        return n/2
    @staticmethod
    def median(n):
        if(n<=0 or n%1!=0):
            raise ValueError("n must be a positive integer")
        return n/2
    @staticmethod
    def mode(n):
        if(n<=0 or n%1!=0):
            raise ValueError("n must be a positive integer")
        if(n==1):
            return 1
        return n/2
    @staticmethod
    def variance(n):
        if(n<=0 or n%1!=0):
            raise ValueError("n must be a positive integer")
        return n/12
    @staticmethod
    def stddev(n):
        if(n<=0 or n%1!=0):
            raise ValueError("n must be a positive integer")
        return math.sqrt(n/12)
    @staticmethod
    def skewness(n):
        if(n<=0 or n%1!=0):
            raise ValueError("n must be a positive integer")
        return 0
