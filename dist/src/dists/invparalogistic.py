'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.genbetaprime.genbetaprime as genbetaprime

class invparalogistic(Distribution): #inverse paralogistic
    @staticmethod
    def random(bb):
        return genbetaprime.random(0,1,bb,1,bb)
