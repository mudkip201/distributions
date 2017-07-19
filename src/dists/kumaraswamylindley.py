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

class kumaraswamylindley(Distribution):
    @staticmethod
    def pdf(a,b,t,x):
        return a*b*t**2/(t+1)*(1+x)*math.exp(-t*x)*math.pow(1-math.exp(-t*x)*(1+(t*x/(t+1))),a-1)*math.pow(1-math.pow(1-math.exp(-t*x)*(1+t*x/(t+1)),a),b-1)
    @staticmethod
    def cdf(a,b,t,x):
        return 1-math.pow(1-math.pow(1-math.exp(-t*x)*(1+t*x/(t+1)),a),b)
