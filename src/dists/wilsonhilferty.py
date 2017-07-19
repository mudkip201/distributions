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
import dists.amoroso.amoroso as amoroso

class wilsonhilferty(Distribution):
    @staticmethod
    def random(theta,aa):
        return amoroso.random(0,theta,aa,3)
    @staticmethod
    def pdf(theta,aa,x):
        return amoroso.pdf(0,theta,aa,3,x)
    @staticmethod
    def mean(theta,aa):
        return amoroso.mean(0,theta,aa,3)
    @staticmethod
    def mode(theta,aa):
        return amoroso.mode(0,theta,aa,3)
    @staticmethod
    def variance(theta,aa):
        return amoroso.variance(0,theta,aa,3)
    @staticmethod
    def stddev(theta,aa):
        return amoroso.stddev(0,theta,aa,3)