'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.power.power as power

class wedge(Distribution):
    @staticmethod
    def random(a,s):
        return power.random(a,s,2)