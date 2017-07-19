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

class bates(Distribution):
    @staticmethod
    def random(n):
        if(n%1!=0 or n<1):
            raise ValueError("n must be a positive integer greater than 0")
        avg=0
        for _ in range(n):
            avg+=r.random()
        return avg/n
    @staticmethod
    def kurtosis(n):
        if(n%1!=0 or n<1):
            raise ValueError("n must be a positive integer greater than 0")
        return -6/(5*n)
    @staticmethod
    def mean(n):
        if(n%1!=0 or n<1):
            raise ValueError("n must be a positive integer greater than 0")
        return 1/2
    @staticmethod
    def variance(n):
        if(n%1!=0 or n<1):
            raise ValueError("n must be a positive integer greater than 0")
        return 1/(12*n)
    @staticmethod
    def stddev(n):
        if(n%1!=0 or n<1):
            raise ValueError("n must be a positive integer greater than 0")
        return math.sqrt(1/(12*n))
    @staticmethod
    def skewness(n):
        if(n%1!=0 or n<1):
            raise ValueError("n must be a positive integer greater than 0")
        return 0