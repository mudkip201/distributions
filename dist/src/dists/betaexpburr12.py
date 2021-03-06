'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.beta.beta as beta

class betaexpburr12(Distribution): #beta-exponential burr XII
    @staticmethod
    def random(a,b,c,d,k):
        u=beta.random(a,b)
        return math.pow(math.pow(1-math.pow(u,1/b),-1/k)-1,1/c)