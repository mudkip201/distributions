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

class gengamma(Distribution): #generalized gamma
    @staticmethod
    def random(a,c):
        return st.gengamma.rvs(a,c)