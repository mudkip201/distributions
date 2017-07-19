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
import dists.normal.normal as normal
import dists.chi2.chi2 as chi2

class noncentralt(Distribution):
    @staticmethod
    def random(mu,nu):
        return (normal.random(0,1)+mu)/math.sqrt(chi2.random(nu)/nu)