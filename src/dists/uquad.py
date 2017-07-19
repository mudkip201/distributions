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

class uquad(Distribution):
    @staticmethod
    def random(a,b):
        c=12/(math.pow((b-a),3))
        d=(a+b)/2
        return 3*r.random()/c-math.pow(math.pow((d-c),3),1/3.0)+d
    @staticmethod
    def median(a,b):
        c=12/(math.pow((b-a),3))
        d=(a+b)/2
        return 3*(1/2)/c-math.pow(math.pow((d-c),3),1/3.0)+d
    @staticmethod
    def ppf(a,b,q):
        c=12/(math.pow((b-a),3))
        d=(a+b)/2
        return 3*q/c-math.pow(math.pow((d-c),3),1/3.0)+d
