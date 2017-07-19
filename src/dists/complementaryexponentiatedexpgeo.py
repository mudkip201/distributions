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

class complementaryexponentiatedexpgeo(Distribution):
    @staticmethod
    def random(aa,lmbda,theta):
        n=ds.rg0()
        return -math.log(1-math.pow(n/(theta*(1-n)+n),1/aa))/lmbda
    @staticmethod
    def pdf(aa,lmbda,theta,x):
        return aa*lmbda*theta*math.exp(-lmbda*x)*math.pow(1-math.exp(-lmbda*x),aa-1)/math.pow(1-(1-theta)*math.pow(1-math.exp(-lmbda*x),aa),2)
    @staticmethod
    def cdf(aa,lmbda,theta,x):
        return 1-(1-math.pow(1-math.exp(-lmbda*x),aa))/(1-(1-theta)*math.pow(1-math.exp(-lmbda*x),aa))