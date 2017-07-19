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

class wakeby(Distribution):
    @staticmethod
    def random(a,b,c,d,m):
        n=r.random()
        #return -a*pow(1-n,b)+c*pow(1-n,-d)+e
        return m+a*(1-math.pow(1-n,b))-c*(1-math.pow(1-n,-d))
    @staticmethod
    def median(a,b,c,d,m):
        return m+a*(1-math.pow(1/2,b))-c*(1-math.pow(1/2,-d))
    @staticmethod
    def ppf(a,b,c,d,m,q):
        return m+a*(1-math.pow(1-q,b))-c*(1-math.pow(1-q,-d))