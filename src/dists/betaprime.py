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
import dists.gamma.gamma as gamma

class betaprime(Distribution):
    @staticmethod
    def pdf(aa,bb,x):
        return math.pow(x,aa-1)*math.pow(1+x,-aa-bb)/sp.beta(aa,bb)
    @staticmethod
    def cdf(aa,bb,x):
        return sp.betainc(aa,bb,x/(1+x))*sp.beta(aa,bb)
    @staticmethod
    def random(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be greater than 0")
        return gamma.random(aa,1)/gamma.random(bb,1)
    @staticmethod
    def mean(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be greater than 0")
        if(bb>1):
            return aa/(bb-1)
        return None
    @staticmethod
    def mode(aa,bb):
        if(aa<=0 or bb<=0):
            raise ValueError("aa and bb must be greater than 0")
        if(aa>=1):
            return (aa-1)/(bb+1)
        return 0
    @staticmethod
    def variance(aa,bb):
        if(bb>2):
            return aa*(aa+bb-1)/((bb-2)**2*(bb-1)**2)
        return None
    @staticmethod
    def stddev(aa,bb):
        if(bb>2):
            return math.sqrt(betaprime.variance(aa,bb))
        return None
    @staticmethod
    def skewness(aa,bb):
        if(bb>3):
            return 2*(2*aa+bb-1)/(bb-3)*math.sqrt((bb-2)/(aa*(aa+bb-1)))
        return None