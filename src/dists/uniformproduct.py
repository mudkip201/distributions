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

class uniformproduct(Distribution):
    @staticmethod
    def random(n):
        prod=1
        for _ in range(n):
            prod*=r.random()
        return prod