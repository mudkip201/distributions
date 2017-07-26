'''
Created on Jul 26, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class biparabolic(Distribution):
    @staticmethod
    def pdf(a,b,m,x):
        if(a<=x and x<=m):
            return -3/2*1/(m-a)**2/(b-a)*(x**2-2*m*x+(2*m-a)*a)
        if(m<=x and x<=b):
            return -3/2*1/(m-b)**2/(b-a)*(x**2-2*m*x+(2*m-b)*b)
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,b,m):
        return (2*m+3*b+3*a)/8
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a,b,m):
        return (12*(m-a)**2-12*(m-a)*(b-a)+19*(b-a)**2)/320
    @staticmethod
    def stddev(a,b,m):
        return math.sqrt((12*(m-a)**2-12*(m-a)*(b-a)+19*(b-a)**2)/320)
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