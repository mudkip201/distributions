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

class bradford(Distribution):
    @staticmethod
    def random(min_,max_,theta):
        n=r.random()
        return (-min_*math.pow(theta+1,n)+theta*min_+max_*math.pow(theta+1,n)+min_-max_)/theta