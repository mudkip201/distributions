'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
from numpy import random as r
import dists.halfgenpearson7.halfgenpearson7 as halfgenpearson7

class genpearson7(Distribution): #generalized pearson VII
    @staticmethod
    def random(a,s,m,b):
        n=r.random()
        if(n<0.5):
            return -halfgenpearson7.random(a,s,m,b)
        return halfgenpearson7.random(a,s,m,b)
