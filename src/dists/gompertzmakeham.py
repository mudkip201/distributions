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

class gompertzmakeham(Distribution):
    @staticmethod
    def random(lmbda,xi):
        q=ds.rg0()
        return math.log(1-math.log(1-q)/xi)/lmbda
    @staticmethod
    def median(lmbda,xi):
        return math.log(1-math.log(1/2)/xi)/lmbda
    @staticmethod
    def ppf (lmbda,xi,q):
        return math.log(1-math.log(1-q)/xi)/lmbda
