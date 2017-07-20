'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class genlambda(Distribution):
    @staticmethod
    def random(l1,l2,l3,l4):
        n=r.random()
        return l1+(math.pow(n,l3)-math.pow((1-n),l4))/l2
    @staticmethod
    def median(l1,l2,l3,l4):
        return l1+(math.pow(1/2,l3)-math.pow(1/2,l4))/l2
    @staticmethod
    def ppf(l1,l2,l3,l4,q):
        return l1+(math.pow(q,l3)-math.pow((1-q),l4))/l2
