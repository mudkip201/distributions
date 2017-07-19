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

class gandh(Distribution):
    @staticmethod
    def random(g,h):
        n=ds.rg0()
        return math.exp(g*st.norm.ppf(n)-1)*(math.exp((h*math.pow(st.norm.ppf(n),2))/2)/g)
    @staticmethod
    def median(g,h):
        return math.exp(g*st.norm.ppf(1/2)-1)*(math.exp((h*math.pow(st.norm.ppf(1/2),2))/2)/g)
    @staticmethod
    def ppf(g,h,q):
        return math.exp(g*st.norm.ppf(q)-1)*(math.exp((h*math.pow(st.norm.ppf(q),2))/2)/g)