'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Dist
import math
import dists.beta.beta as beta


class baldingnichols(Dist):
    @staticmethod
    def random(f,p):
        if(f<=0 or f>=1 or p<=0 or p>=1):
            raise ValueError("all inputs must be between 0 and 1 exclusive")
        return beta.random(p*(1-f)/f,(1-p)*(1-f)/f)
    @staticmethod
    def mean(f,p):
        if(f<=0 or f>=1 or p<=0 or p>=1):
            raise ValueError("all inputs must be between 0 and 1 exclusive")
        return p
    @staticmethod
    def mode(f,p):
        if(f<=0 or f>=1 or p<=0 or p>=1):
            raise ValueError("all inputs must be between 0 and 1 exclusive")
        return (f-(1-f)*p)/(3*f-1)
    @staticmethod
    def variance(f,p):
        if(f<=0 or f>=1 or p<=0 or p>=1):
            raise ValueError("all inputs must be between 0 and 1 exclusive")
        return f*p*(1-p)
    @staticmethod
    def stddev(f,p):
        if(f<=0 or f>=1 or p<=0 or p>=1):
            raise ValueError("all inputs must be between 0 and 1 exclusive")
        return math.sqrt(f*p*(1-p))
    @staticmethod
    def skewness(f,p):
        if(f<=0 or f>=1 or p<=0 or p>=1):
            raise ValueError("all inputs must be between 0 and 1 exclusive")
        return 2*f*(1-2*p)/((1+f)*math.sqrt(f*(1-p)*p))