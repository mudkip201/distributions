'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
from numpy import random as r

class uniformproduct(Distribution):
    @staticmethod
    def random(n):
        prod=1
        for _ in range(n):
            prod*=r.random()
        return prod