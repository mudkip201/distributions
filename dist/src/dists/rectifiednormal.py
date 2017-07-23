'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.normal.normal as normal

class rectifiednormal(Distribution):
    @staticmethod
    def random(mu,sigma):
        n=normal.random(mu,sigma)
        if(n<0):
            return 0
        return n