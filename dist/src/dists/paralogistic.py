'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.genbetaprime.genbetaprime as genbetaprime

class paralogistic(Distribution):
    @staticmethod
    def random(bb):
        return genbetaprime.random(0,1,1,bb,bb)