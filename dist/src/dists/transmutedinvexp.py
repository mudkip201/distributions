'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class transmutedinvexp(Distribution): #transmuted inverse exponential
    @staticmethod
    def random(a,l):
        u=ds.rg0()
        i=2*l/((1+l)-math.sqrt((1+l)**2-4*l*u))
        return a/math.log(i)