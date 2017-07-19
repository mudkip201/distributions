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

class expgeninvexp(Distribution): #exponentiated generalized inverse exponential
    @staticmethod
    def random(a,b,l):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(ds.rg0(),1/b),1/a)),-1)
    @staticmethod
    def median(a,b,l):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(1/2,1/b),1/a)),-1)
    @staticmethod
    def ppf(a,b,l,q):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(q,1/b),1/a)),-1)