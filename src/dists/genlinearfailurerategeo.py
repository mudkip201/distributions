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

class genlinearfailurerategeo(Distribution): #generalized linear failure rate-geometric
    @staticmethod
    def random(a,b,aa,p):
        u=ds.rg0()
        return -a/b+math.sqrt(a**2-2*b*math.log(1-math.pow((1-p)*u)/(1-p*u),1/aa))/b
    @staticmethod
    def median(a,b,aa,p):
        return -a/b+math.sqrt(a**2-2*b*math.log(1-math.pow((1-p)/2)/(1-p/2),1/aa))/b
    @staticmethod
    def ppf(a,b,aa,p,q):
        return -a/b+math.sqrt(a**2-2*b*math.log(1-math.pow((1-p)*q)/(1-p*q),1/aa))/b