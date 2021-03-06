'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.amoroso.amoroso as amoroso

class gammaexp(Distribution): #gamma-exponential
    @staticmethod
    def random(nu,lmbda,a):
        return -math.log(amoroso.random(0,math.exp(nu),a,1/lmbda))