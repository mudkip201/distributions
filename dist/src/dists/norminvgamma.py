'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.invgamma.invgamma as invgamma
import dists.normal.normal as normal

class norminvgamma(Distribution): #normal-inverse gamma
    @staticmethod
    def random(aa,bb,mu,lmbda):
        n=invgamma.random(aa,bb)
        return normal.random(mu,math.sqrt(n/lmbda))
