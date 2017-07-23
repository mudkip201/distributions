'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class skewlaplace(Distribution):
    @staticmethod
    def random(a,b,c):
        n=r.random()
        if(n<a):
            return a-b*math.log(b/(n*(b+c)))
        return a-b*math.log(-((n-1)*(b+c))/b)
    @staticmethod
    def median(a,b,c):
        if(1/2<a):
            return a-b*math.log(b/((b+c)/2))
        return a-b*math.log((b+c)/(2*b))
    @staticmethod
    def ppf(a,b,c,q):
        if(q<a):
            return a-b*math.log(b/(q*(b+c)))
            return a-b*math.log(-((q-1)*(b+c))/b)