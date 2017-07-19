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

class kumaraswamyburr12(Distribution):
    @staticmethod
    def pdf(a,b,c,k,s,x):
        return a*b*c*k*math.pow(s,-c)*math.pow(x,c-1)*math.pow(1+math.pow(x/s,c),-k-1)*math.pow(1-math.pow(1+math.pow(x/s,c),-k),a-1)*math.pow(1-math.pow(1-math.pow(1+math.pow(x/s,c),-k),a),b-1)
    @staticmethod
    def cdf(a,b,c,k,s,x):
        return 1-math.pow(1-math.pow(1-math.pow(1+math.pow(x/s,c),-k),a),b)
    @staticmethod
    def random(a,b,c,k,s):
        u=ds.rg0()
        return s*math.pow(math.pow(1-math.pow(1-math.pow(1-u,1/b),1/a),-1/k)-1,1/c)
    @staticmethod
    def median(a,b,c,k,s):
        return s*math.pow(math.pow(1-math.pow(1-math.pow(1/2,1/b),1/a),-1/k)-1,1/c)
    @staticmethod
    def ppf(a,b,c,k,s,q):
        return s*math.pow(math.pow(1-math.pow(1-math.pow(1-q,1/b),1/a),-1/k)-1,1/c)