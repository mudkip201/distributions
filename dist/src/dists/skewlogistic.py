'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r


class skewlogistic(Distribution):
    @staticmethod
    def random(a,b,c):
        n=r.random()
        return math.log(math.pow(1-n,-1/c))*(a*math.pow(math.log(1-n),1/c)-b)
    @staticmethod
    def median(a,b,c):
        return math.log(math.pow(-1/2,-1/c))*(a*math.pow(math.log(1/2),1/c)-b)
    @staticmethod
    def ppf(a,b,c,q):
        return math.log(math.pow(1-q,-1/c))*(a*math.pow(math.log(1-q),1/c)-b)