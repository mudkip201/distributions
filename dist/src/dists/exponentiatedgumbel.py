'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class exponentiatedgumbel(Distribution):
    @staticmethod
    def random(a,m,s):
        return -s*math.log(-math.log(1-math.pow(1-ds.rg0(),1/a)))+m
    @staticmethod
    def cdf(a,m,s,x):
        return 1-math.pow(1-math.exp(-math.exp(-(x-m)/s)),a)
    @staticmethod
    def median(a,m,s):
        return -s*math.log(-math.log(1-math.pow(1/2,1/a)))+m
    @staticmethod
    def ppf(a,m,s,q):
        return -s*math.log(-math.log(1-math.pow(1-q,1/a)))+m