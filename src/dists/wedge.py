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
import dists.power.power as power

class wedge(Distribution):
    @staticmethod
    def random(a,s):
        return power.random(a,s,2)