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
import dists.normal.normal as normal

class truncnormr(Distribution): #right-truncated normal #truncated on the right
    @staticmethod
    def random(b,mu,sigma):
        n=normal.random(mu,sigma**2)
        while n>b:
            n=normal.random(mu,sigma**2)
        return n