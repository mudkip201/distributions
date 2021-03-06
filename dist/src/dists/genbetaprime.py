'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.gamma.gamma as gamma

class genbetaprime(Distribution): #generalized beta-prime
    @staticmethod
    def random(a,s,aa,gmma,bb):
        return a+s*math.pow((gamma.random(aa,1)/gamma.random(gmma,1)),1/bb)
