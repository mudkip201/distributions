'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.normal.normal as normal
import dists.cauchy.cauchy as cauchy

class voigt(Distribution):
    @staticmethod
    def random(a,s,sigma):
        return normal.random(0,sigma)+cauchy.random(a,s)
