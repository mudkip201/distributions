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

class exponentiatedgompertz(Distribution):
    @staticmethod
    def random(a,l,t):
        return math.log(1-math.log(1-math.pow(ds.rg0(),1/t))/l)/a
    @staticmethod
    def median(a,l,t):
        return math.log(1-math.log(1-math.pow(1/2,1/t))/l)/a
    @staticmethod
    def ppf(a,l,t,q):
        return math.log(1-math.log(1-math.pow(q,1/t))/l)/a