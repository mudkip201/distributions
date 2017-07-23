'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class transmutedgenmodweibull(Distribution): #transmuted generalized modified weibull
    @staticmethod
    def pdf(aa,gmma,eta,lmbda,phi,x):
        u=aa*phi*math.pow(x,gmma-1)*(x*lmbda+gmma)
        u*=math.exp(x*lmbda-aa*math.pow(x,gmma)*math.exp(x,lmbda))
        u*=math.pow(1-math.exp(-aa*math.pow(x,gmma)*math.exp(x*lmbda)),phi-1)
        return u*(1+eta-2*eta*math.pow(1-math.exp(-aa*math.pow(x,gmma)*math.exp(x*lmbda)),phi))
    @staticmethod
    def cdf(aa,gmma,eta,lmbda,phi,x):
        u=math.pow(1-math.exp(-aa*math.pow(x,gmma)*math.exp(x*lmbda)),phi)*(1+eta-eta*math.pow(1-math.exp(-aa*math.pow(x,gmma)*math.exp(x*lmbda)),phi))