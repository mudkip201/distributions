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


class benford(Distribution):
    @staticmethod
    def random(b):
        return 1/(math.pow(b,1-r.random()))