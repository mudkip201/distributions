'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.prentice.prentice as prentice

class reversedburr2(Distribution):
    @staticmethod
    def random(g):
        return prentice.random(0,-1,1,g)