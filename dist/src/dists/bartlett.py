'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Dist
import dists.poisson.poisson as poisson
import dists.geometric.geometric as geometric

class bartlett(Dist):
    @staticmethod
    def random(a,q):
        return poisson.random(a)*geometric.random(q)