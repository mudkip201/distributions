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

class gentypeIextremevalue(Distribution): #generalized type I extreme value
    @staticmethod
    def random(a,b,m,s):
        return -s*math.log(b/s*(math.pow(ds.rg0(),-1/a)-1))+m
    @staticmethod
    def median(a,b,m,s):
        return -s*math.log(b/s*(math.pow(1/2,-1/a)-1))+m
    @staticmethod
    def ppf(a,b,m,s,q):
        return -s*math.log(b/s*(math.pow(q,-1/a)-1))+m
