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
import dists.genpearson7.genpearson7 as genpearson7

class meridian(Distribution):
    @staticmethod
    def random(a,s):
        return genpearson7.random(a,s,2,1)