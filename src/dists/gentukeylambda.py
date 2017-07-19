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

class gentukeylambda(Distribution): #generalized tukey-lambda
    @staticmethod
    def random(a,d):
        u=r.random()
        return (math.pow(u,a-d)-1)/(a-d)-(math.pow(1-u,a+d)-1)/(a+d)
    @staticmethod
    def median(a,d):
        return (math.pow(1/2,a-d)-1)/(a-d)-(math.pow(1/2,a+d)-1)/(a+d)
    @staticmethod
    def ppf(a,d,q):
        return (math.pow(q,a-d)-1)/(a-d)-(math.pow(1-q,a+d)-1)/(a+d)
