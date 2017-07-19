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

class extremevaluemin(Distribution):
    @staticmethod
    def random(a,b):
        return -a+b*math.log(math.log(1/ds.rg0()))
