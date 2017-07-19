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

class anglit(Distribution):
    @staticmethod
    def pdf(x):
        return math.cos(2*x)
    @staticmethod
    def cdf(x):
        return math.sin(x+math.pi/4)**2
    @staticmethod
    def random():
        n=r.random()
        return math.asin(math.sqrt(n))-math.pi/4
    @staticmethod
    def median():
        return math.asin(math.sqrt(1/2))-math.pi/4
    @staticmethod
    def ppf(q):
        return math.asin(math.sqrt(q))-math.pi/4