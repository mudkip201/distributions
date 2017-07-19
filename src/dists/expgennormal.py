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

class expgennormal(Distribution): #exponentiated generalized normal
    @staticmethod
    def random(a,b,m,s):
        u=ds.rg0()
        return st.norm.ppf(s*(1-math.pow(1-math.pow(u,1/b),1/a))+m,0,1)
    @staticmethod
    def median(a,b,m,s):
        return st.norm.ppf(s*(1-math.pow(1-math.pow(1/2,1/b),1/a))+m,0,1)
    @staticmethod
    def ppf(a,b,m,s,q):
        return st.norm.ppf(s*(1-math.pow(1-math.pow(q,1/b),1/a))+m,0,1)