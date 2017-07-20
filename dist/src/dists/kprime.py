'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.normal.normal as normal
import dists.chi2.chi2 as chi2

class kprime(Distribution):
    @staticmethod
    def random(a,b,nu1,nu2):
        return (b*normal.random(0,1)+a*math.sqrt(chi2.random(nu1)/nu1))/math.sqrt(chi2.random(nu2)/nu2)