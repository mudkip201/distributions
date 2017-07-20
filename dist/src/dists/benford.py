'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r


class benford(Distribution):
    @staticmethod
    def random(b):
        return 1/(math.pow(b,1-r.random()))