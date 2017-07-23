'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class slope(Distribution):
    @staticmethod
    def random(aa):
        if aa==1:
            return 1
        return (-aa+math.sqrt(math.pow(aa,2)+4*r.random()*(1-aa)))/(2*(1-aa))
    @staticmethod
    def median(aa):
        if aa==1:
            return 1
        return (-aa+math.sqrt(math.pow(aa,2)+2*(1-aa)))/(2*(1-aa))
    @staticmethod
    def ppf(aa,q):
        if aa==1:
            return 1
        return (-aa+math.sqrt(math.pow(aa,2)+4*q*(1-aa)))/(2*(1-aa))