'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.stdbetaprime.stdbetaprime as stdbetaprime

class stdprentice(Distribution): #standard prentice
    @staticmethod
    def random(a,g):
        return math.log(stdbetaprime.random(a,g))