'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.genbetaprime.genbetaprime as genbetaprime

class moffat(Distribution):
    @staticmethod
    def random(s,g):
        return genbetaprime.random(0,s,1,g,2)