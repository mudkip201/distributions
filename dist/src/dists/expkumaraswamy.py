'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class expkumaraswamy(Distribution): #exponentiated kumaraswamy
    @staticmethod
    def random(a,b,g):
        return math.pow(1-math.pow(1-math.pow(ds.rg0(),1/g),1/b),1/a)
    @staticmethod
    def median(a,b,g):
        return math.pow(1-math.pow(1-math.pow(1/2,1/g),1/b),1/a)
    @staticmethod
    def ppf(a,b,g,q):
        return math.pow(1-math.pow(1-math.pow(q,1/g),1/b),1/a)
