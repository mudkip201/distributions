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
import dists.gammaexp.gammaexp as gammaexp

class bramwellholdsworthpinton(Distribution):
    @staticmethod
    def random(nu,lmbda):
        return gammaexp.random(nu,lmbda,math.pi/2)