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

class kumaraswamy4(Distribution):
    @staticmethod
    def random(a,b,c,d):
        '''alpha, beta, min, max'''
        n=ds.rg0()
        while n==1:
            n=ds.rg0()
        return c+(d-c)*math.pow(1-math.pow(n,1/b),1/a)
    @staticmethod
    def median(a,b,c,d):
        '''alpha, beta, min, max'''
        return c+(d-c)*math.pow(1-math.pow(1/2,1/b),1/a)
    @staticmethod
    def ppf(a,b,c,d,q):
        '''alpha, beta, min, max'''
        return c+(d-c)*math.pow(1-math.pow(q,1/b),1/a)