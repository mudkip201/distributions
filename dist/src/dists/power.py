'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class power(Distribution):
    @staticmethod
    def random(c,a,b):
        n=ds.rg0()
        return a+b*math.pow(n,1/c)
    @staticmethod
    def median(c,a,b):
        return a+b*math.pow(1/2,1/c)
    @staticmethod
    def ppf(c,a,b,q):
        return a+b*math.pow(q,1/c)
