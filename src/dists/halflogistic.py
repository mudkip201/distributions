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
import dists.logistic.logistic as logistic

class halflogistic(Distribution):
    @staticmethod
    def random(mu,s):
        n=logistic.random(mu,s)
        while(n<0):
            n=logistic.random(mu,s)
        return n
