'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.gamma.gamma as gamma

class wein(Distribution):
    @staticmethod
    def random(t):
        return gamma.random(t,4)