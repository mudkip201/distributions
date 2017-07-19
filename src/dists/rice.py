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
import dists.chi2.chi2 as chi2
import dists.poisson.poisson as poisson

class rice(Distribution):
    @staticmethod
    def random(nu,sigma):
        return sigma*math.sqrt(chi2.random(2*poisson.random(nu**2/(2*sigma**2))+2))