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

class stacy(Distribution):
    @staticmethod
    def random(theta,aa,bb):
        return amoroso.random(0,theta,aa,bb)
    @staticmethod
    def pdf(theta,aa,bb,x):
        return amoroso.pdf(0,theta,aa,bb,x)
    @staticmethod
    def mean(theta,aa,bb):
        return amoroso.mean(0,theta,aa,bb)
    @staticmethod
    def mode(theta,aa,bb):
        return amoroso.mode(0,theta,aa,bb)
    @staticmethod
    def variance(theta,aa,bb):
        return amoroso.variance(0,theta,aa,bb)
    @staticmethod
    def stddev(theta,aa,bb):
        return amoroso.stddev(0,theta,aa,bb)