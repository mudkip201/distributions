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

class transmutedexponentiatedexp(Distribution): #transmuted exponentiated exponential
    @staticmethod
    def random(a,g,l):
        i=((1+l)-math.sqrt((1+l)**2-4*l*ds.rg0()))/(2*l)
        return -math.log(1-math.pow(i,1/a))/g
    @staticmethod
    def median(a,g,l):
        i=((1+l)-math.sqrt((1+l)**2-2*l))/(2*l)
        return -math.log(1-math.pow(i,1/a))/g
    @staticmethod
    def ppf(a,g,l,q):
        i=((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return -math.log(1-math.pow(i,1/a))/g
