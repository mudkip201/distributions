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

class weibullweibull(Distribution):
    @staticmethod
    def random(a,b,g,l):
        return math.pow(math.log(math.pow(-math.log(1-ds.rg0())/a,1/b)+1)/l,1/g)
    @staticmethod
    def median(a,b,g,l):
        return math.pow(math.log(math.pow(-math.log(1/2)/a,1/b)+1)/l,1/g)
    @staticmethod
    def ppf(a,b,g,l,q):
        return math.pow(math.log(math.pow(-math.log(1-q)/a,1/b)+1)/l,1/g)