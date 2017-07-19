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

class arcsin(Distribution):
    @staticmethod
    def pdf(x):
        if(x<0 or x>1):
            raise ValueError("Must be between 0 and 1")
        return 1/(math.pi*math.sqrt(x*(1-x)))
    @staticmethod
    def cdf(x):
        if(x<0 or x>1):
            raise ValueError("Must be between 0 and 1")
        return math.asin(2*x-1)/math.pi+1/2
    @staticmethod
    def random():
        return math.pow(math.sin(math.pi/2*r.random()),2)
    @staticmethod
    def mean():
        return 1/2
    @staticmethod
    def median():
        return 1/2
    @staticmethod
    def variance():
        return 1/8
    @staticmethod
    def stddev():
        return math.sqrt(1/8)
    @staticmethod
    def ppf(q):
        math.pow(math.sin(math.pi/2*q),2)