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

class skewgent(Distribution):
    @staticmethod
    def pdf(k,n,l,s2,x):
        S=math.sqrt(1+3*l**2-4l**2*math.pow(sp.beta(2/k,(n-1)/k),2)**2/sp.beta(1/k,n/k)/sp.beta(3/k,(n-2)/k))
        C=1/2*k*math.pow(sp.beta(1/k,n/k),-3/2)*math.sqrt(sp.beta(3/k,(n-2)/k))*S/math.sqrt(s2)
        theta=math.pow(k/(n-2),1/k)*math.sqrt(sp.beta(1/k,n/k))/math.sqrt(sp.beta(3/k,(n-2)/k))/S
        if(x>=0):
            return C*math.pow(1+(k/(n-2))*math.pow(theta,-k)*math.pow(1+l,-k)*math.pow(abs(x/math.sqrt(s2)),k),-(n+1)/k)
        return C*math.pow(1+(k/(n-2))*math.pow(theta,-k)*math.pow(1-l,-k)*math.pow(abs(x/math.sqrt(s2)),k),-(n+1)/k)
    @staticmethod
    def mean(k,n,l,s2):
        S=math.sqrt(1+3*l**2-4l**2*math.pow(sp.beta(2/k,(n-1)/k),2)**2/sp.beta(1/k,n/k)/sp.beta(3/k,(n-2)/k))
        return 2*l/S*sp.beta(2/k,(n-1)/k)/math.sqrt(sp.beta(1/k,n/k))/math.sqrt(3/k,(n-2)/k)*math.sqrt(s2)
    @staticmethod
    def variance(k,n,l,s2):
        return s2
    @staticmethod
    def stddev(k,n,l,s2):
        return math.sqrt(s2)
