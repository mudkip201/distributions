'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.gammaexp.gammaexp as gammaexp

class bramwellholdsworthpinton(Distribution):
    @staticmethod
    def random(nu,lmbda):
        return gammaexp.random(nu,lmbda,math.pi/2)