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

class twosidedpower(Distribution):
    @staticmethod
    def random(d,theta):
        n=ds.rg0()
        if(n<=theta):
            return theta*math.pow(n/theta,1/d)
        return 1-(1-theta)*math.pow((1-n)/(1-d),1/d)
    @staticmethod
    def median(d,theta):
        if(1/2<=theta):
            return theta*math.pow(1/(2*theta),1/d)
        return 1-(1-theta)*math.pow(1/(2*(1-d)),1/d)
    @staticmethod
    def ppf(d,theta,q):
        if(q<=theta):
            return theta*math.pow(q/theta,1/d)
        return 1-(1-theta)*math.pow((1-q)/(1-d),1/d)