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

class cosine(Distribution):
    @staticmethod
    def random():
        return st.cosine.rvs()
    @staticmethod
    def pdf(x):
        return (1+math.cos(x))/(2*math.pi)
    @staticmethod
    def cdf(x):
        return (math.pi+x+math.sin(x))/(2*math.pi)