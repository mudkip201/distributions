'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.rayleigh.rayleigh as rayleigh
import dists.lognormal.lognormal as lognormal

class suzuki(Distribution):
    @staticmethod
    def random(mu,nu):
        return rayleigh.random(lognormal.random(0,mu,nu))
    @staticmethod
    def kurtosis(mu,nu):
        return (32*math.exp(6*nu**2)+24*math.exp(nu**2)*math.pi-24*math.exp(3*nu**2)*math.pi-3*math.pi**2)/((4*math.exp(nu**2)-math.pi)**2)
    @staticmethod
    def mean(mu,nu):
        return math.exp(mu+nu**2/2)*math.sqrt(math.pi/2)
    @staticmethod
    def variance(mu,nu):
        return math.exp(2*mu+nu**2)*(2*math.exp(nu**2)-math.pi/2)
    @staticmethod
    def stddev(mu,nu):
        return math.sqrt(math.exp(2*mu+nu**2)*(2*math.exp(nu**2)-math.pi/2))
    @staticmethod
    def skewness(mu,nu):
        return 2*math.sqrt(math.pi)*(-6*math.exp(nu**2)+3*math.exp(3*nu**2)+math.pi)/math.pow(4*math.exp(nu**2)-math.pi,3/2)
