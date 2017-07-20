'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.gammaexp.gammaexp as gammaexp

class gengumbel(Distribution): #generalized gumbel
    @staticmethod
    def random(u,lmbda,n):
        return gammaexp.random(u-lmbda*math.log(n),lmbda,n)