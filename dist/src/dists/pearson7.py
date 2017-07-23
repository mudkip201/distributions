'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.normal.normal as normal
import dists.gamma.gamma as gamma

class pearson7(Distribution):
    @staticmethod
    def random(s,m):
        return normal.random(0,s)/math.sqrt(gamma.random(1/2,m-1/2))
