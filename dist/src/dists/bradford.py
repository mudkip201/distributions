'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class bradford(Distribution):
    @staticmethod
    def random(min_,max_,theta):
        n=r.random()
        return (-min_*math.pow(theta+1,n)+theta*min_+max_*math.pow(theta+1,n)+min_-max_)/theta