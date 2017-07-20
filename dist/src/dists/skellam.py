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
import dists.poisson.poisson as poisson

class skellam(Distribution):
    @staticmethod
    def random(mu1,mu2):
        return poisson.random(mu1)-poisson.random(mu2)
    @staticmethod
    def kurtosis(mu1,mu2):
        return 1/(mu1+mu2)
    @staticmethod
    def mean(mu1,mu2):
        return mu1-mu2
    @staticmethod
    def variance(mu1,mu2):
        return mu1+mu2
    @staticmethod
    def stddev(mu1,mu2):
        return math.sqrt(mu1+mu2)
    @staticmethod
    def skewness(mu1,mu2):
        return (mu1-mu2)/math.pow(mu1+mu2,3/2)
