'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.gamma.gamma as gamma
import dists.normal.normal as normal

class normalgamma(Distribution):
    @staticmethod
    def random(aa,bb,mu,lmbda):
        n=gamma.random(aa,bb)
        return normal.random(mu,math.sqrt(1/(lmbda*n)))