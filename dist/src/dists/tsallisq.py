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

class tsallisq(Distribution):
    @staticmethod
    def random(lmbda,q1):
        x=ds.rg0()
        if(q1==1):
            return -math.log(1-x)/lmbda
        return -(1-math.pow(1-x,1/(-2+q1))+math.pow(1-x,1/(-2+q1))*x)/((-1+q1)*lmbda)
    @staticmethod
    def median(lmbda,q1):
        if(q1==1):
            return -math.log(1/2)/lmbda
        return -(1-math.pow(1/2,1/(-2+q1))+math.pow(1/2,1/(-2+q1))/2)/((-1+q1)*lmbda)
    @staticmethod
    def ppf(lmbda,q1,q):
        if(q1==1):
            return -math.log(1-q)/lmbda
        return -(1-math.pow(1-q,1/(-2+q1))+math.pow(1-q,1/(-2+q1))*q)/((-1+q1)*lmbda)