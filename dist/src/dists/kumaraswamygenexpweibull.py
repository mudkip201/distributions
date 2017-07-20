'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class kumaraswamygenexpweibull(Distribution): #kumaraswamy generalized exponentiated weibull
    @staticmethod
    def random(a,b,aa,bb):
        return -math.log(1-math.pow(1-math.pow(1-ds.rg0(),1/b),1/a))/(aa+bb)
    @staticmethod
    def median(a,b,aa,bb):
        return -math.log(1-math.pow(1-math.pow(1/2,1/b),1/a))/(aa+bb)
    @staticmethod
    def ppf(a,b,aa,bb,q):
        return -math.log(1-math.pow(1-math.pow(1-q,1/b),1/a))/(aa+bb)