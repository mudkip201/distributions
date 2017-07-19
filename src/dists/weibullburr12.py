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

class weibullburr12(Distribution):
    @staticmethod
    def random(a,b,c,k,s):
        return s*math.pow(math.pow(math.pow((-math.log(1-ds.rg0())/a),1/b)+1,1/k)-1,1/c)
    @staticmethod
    def median(a,b,c,k,s):
        return s*math.pow(math.pow(math.pow((-math.log(1/2)/a),1/b)+1,1/k)-1,1/c)
    @staticmethod
    def ppf(a,b,c,k,s,q):
        return s*math.pow(math.pow(math.pow((-math.log(1-q)/a),1/b)+1,1/k)-1,1/c)