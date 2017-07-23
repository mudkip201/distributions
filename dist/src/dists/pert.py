'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class pert(Distribution): #PERT distribution
    @staticmethod
    def pdf(minn, maxx, c,l, x):
        if(minn<x and x<maxx):
            return math.pow(maxx-minn,-1-l)*math.pow(maxx-x,l*(maxx-c)/(maxx-minn))*math.pow(x-minn,l*(c-minn)/(maxx-minn))/sp.beta(1+l*(c-minn)/(maxx-minn),1+l*(maxx-c)/(maxx-minn))
        return 0
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(minn, maxx, c, l):
        return (maxx+minn+c*l)/(2+l)
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(maxx,minn,c,l):
        return (maxx-minn-c*l+maxx*l)*(maxx+c*l-minn*(1+l))/((2+l)**2*(3+l))
    @staticmethod
    def stddev(maxx,minn,c,l):
        return math.sqrt((maxx-minn-c*l+maxx*l)*(maxx+c*l-minn*(1+l))/((2+l)**2*(3+l)))
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass