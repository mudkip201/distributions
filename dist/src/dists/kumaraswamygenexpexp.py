'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class kumaraswamygenexpexp(Distribution): #kumaraswamy generalized exponentiated exponential
    @staticmethod
    def random(a,b,c,l):
        return 1/(l*math.log(1-math.pow(1-math.pow(1-ds.rg0(),1/b),1/(a*c))))
    @staticmethod
    def median(a,b,c,l):
        return 1/(l*math.log(1-math.pow(1-math.pow(1/2,1/b),1/(a*c))))
    @staticmethod
    def ppf(a,b,c,l,q):
        return 1/(l*math.log(1-math.pow(1-math.pow(1-q,1/b),1/(a*c))))
