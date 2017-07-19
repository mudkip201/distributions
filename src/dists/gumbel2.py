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

class gumbel2(Distribution):
    @staticmethod
    def random(a,b):
        return math.log(ds.rg0())/(a*b)
    @staticmethod
    def median(a,b):
        return math.log(1/2)/(a*b)
    @staticmethod
    def ppf(a,b,q):
        return math.log(q)/(a*b)
