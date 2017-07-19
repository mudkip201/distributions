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

class qgaussian(Distribution):
    @staticmethod
    def random(bb,q1,mu):
        u=ds.rg0()
        return mu+(math.sqrt(-2*ds.qlog((1+q1)/(3-q1),u))*math.cos(2*math.pi*u))/math.sqrt(bb*(3-q1))
    @staticmethod
    def median(bb,q1,mu):
        return mu+(math.sqrt(-2*ds.qlog((1+q1)/(3-q1),1/2))*math.cos(2*math.pi/2))/math.sqrt(bb*(3-q1))
    @staticmethod
    def ppf(bb,q1,mu,q):
        return mu+(math.sqrt(-2*ds.qlog((1+q1)/(3-q1),q))*math.cos(2*math.pi*q))/math.sqrt(bb*(3-q1))