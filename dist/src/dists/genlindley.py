'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import dists.gamma.gamma as gamma

class genlindley(Distribution): #generalized lindley
    @staticmethod
    def random(a,t,b):
        u=ds.rg0()
        v=gamma.random(a,t)
        vv=gamma.random(a+1,t)
        if(u<=(t/(b+t))):
            return v
        return vv