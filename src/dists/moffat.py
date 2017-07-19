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
import dists.genbetaprime.genbetaprime as genbetaprime

class moffat(Distribution):
    @staticmethod
    def random(s,g):
        return genbetaprime.random(0,s,1,g,2)