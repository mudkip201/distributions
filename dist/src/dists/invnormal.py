'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import dists.normal.normal as normal

class invnormal(Distribution): #inverse normal
    @staticmethod
    def random(mu,lmbda):
        nu=normal.random(0,1)
        y=nu*nu
        x=mu+(mu*mu*y/(2*lmbda))-(mu/(2*lmbda))*math.sqrt(4*mu*lmbda*y+mu*mu*y*y)
        z=r.random()
        if(z<=(mu/(mu+x))):
            return x
        return mu*mu/x
