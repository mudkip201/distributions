'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class newmodweibull(Distribution):
    @staticmethod
    def pdf(a,b,l,x):
        return a*(b+l*x)*math.pow(x,b-1)*math.exp(l*x)*math.exp(-a*math.pow(x,b)*math.exp(x*l))
    @staticmethod
    def cdf(a,b,l,x):
        return 1-math.exp(-a*math.pow(x,b)*math.exp(x*l))