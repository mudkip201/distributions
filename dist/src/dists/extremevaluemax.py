'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class extremevaluemax(Distribution):
    @staticmethod
    def random(a,b):
        return a-b*math.log(math.log(1/ds.rg0()))