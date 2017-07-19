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

class pseudoweibull(Distribution):
    @staticmethod
    def random(theta,bb):
        return amoroso.random(0,theta,1+1/bb,bb)
    @staticmethod
    def pdf(theta,bb,x):
        return amoroso.pdf(0,theta,1+1/bb,bb,x)
    @staticmethod
    def mean(theta,bb):
        return amoroso.mean(0,theta,1+1/bb,bb)
    @staticmethod
    def mode(theta,bb):
        return amoroso.mode(0,theta,1+1/bb,bb)
    @staticmethod
    def variance(theta,bb):
        return amoroso.variance(0,theta,1+1/bb,bb)
    @staticmethod
    def stddev(theta,bb):
        return amoroso.stddev(0,theta,1+1/bb,bb)