'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.normal.normal as normal

class truncnorml(Distribution): #left-truncated normal #truncated on the left
    @staticmethod
    def random(a,mu,sigma):
        n=normal.random(mu,sigma**2)
        while n<a:
            n=normal.random(mu,sigma**2)
        return n
