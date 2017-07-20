'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class genmodweibull(Distribution):
    @staticmethod
    def pdf(aa,gmma,lmbda,phi,x):
        if(x>0):
            u=aa*phi*math.pow(x,gmma-1)*(gmma+lmbda*x)+math.exp(lmbda*x-aa*math.pow(x,gmma)*math.exp(lmbda*x))
            return u/math.pow(1-math.exp(-aa*math.pow(x,gmma)*math.exp(lmbda*x)),1-phi)
    @staticmethod
    def cdf(aa,gmma,lmbda,phi,x):
        return math.pow(1-math.exp(-aa*math.pow(x,gmma)*math.exp(lmbda*x)),phi)