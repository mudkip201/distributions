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
import dists.normal.normal as normal

class slash(Distribution):
    @staticmethod
    def random():
        n=ds.rg0()
        while(n==1):
            n=ds.rg0()
        return normal.random(0,1)/n
    @staticmethod
    def pdf(x):
        return (st.norm.pdf(0)-st.norm.pdf(x))/x**2
    @staticmethod
    def cdf(x):
        if(x!=0):
            return st.norm.cdf(x)-(st.norm.pdf(0)-st.norm.pdf(x))/x
        return 1/2
    @staticmethod
    def median():
        return 0
    @staticmethod
    def mode():
        return 0