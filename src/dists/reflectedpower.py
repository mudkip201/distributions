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

class reflectedpower(Distribution):
    @staticmethod
    def random(c):
        return 1-math.pow(1-r.random(),1/c)
    @staticmethod
    def median(c):
        return 1-math.pow(1/2,1/c)
    @staticmethod
    def ppf(c,q):
        return 1-math.pow(1-q,1/c)