'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class kumaraswamylinearexp(Distribution): #kumaraswamy linear exponential
    @staticmethod
    def random(a,b,l,t):
        return (-l+math.sqrt(l**2-(2*t/a)*math.log(1-(1-math.pow(1-ds.rg0(),1/b)))))/t
    @staticmethod
    def median(a,b,l,t):
        return (-l+math.sqrt(l**2-(2*t/a)*math.log(1-(1-math.pow(1/2,1/b)))))/t
    @staticmethod
    def ppf(a,b,l,t,q):
        return (-l+math.sqrt(l**2-(2*t/a)*math.log(1-(1-math.pow(1-q,1/b)))))/t