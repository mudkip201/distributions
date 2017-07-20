'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.amoroso.amoroso as amoroso

class fishertippett(Distribution):
    @staticmethod
    def random(a,theta,bb):
        return amoroso.random(a,theta,1,bb)
    @staticmethod
    def pdf(a,theta,bb,x):
        return amoroso.pdf(a,theta,1,bb,x)
    @staticmethod
    def mean(a,theta,bb):
        return amoroso.mean(a,theta,1,bb)
    @staticmethod
    def mode(a,theta,bb):
        return amoroso.mode(a,theta,1,bb)
    @staticmethod
    def variance(a,theta,bb):
        return amoroso.variance(a,theta,1,bb)
    @staticmethod
    def stddev(a,theta,bb):
        return amoroso.stddev(a,theta,1,bb)