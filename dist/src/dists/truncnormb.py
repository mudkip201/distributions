'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.normal.normal as normal

class truncnormb(Distribution): #truncated normal (truncated on both ends)
    @staticmethod
    def random(a,b,mu,sigma):
        n=normal.random(mu,sigma**2)
        while (n<a or n>b):
            n=normal.random(mu,sigma**2)
        return n