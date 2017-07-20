'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import numpy as np
import dists.normal.normal as normal

class lognormal(Distribution):
    @staticmethod
    def random(mu,sigma):
        return math.exp(normal.random(mu,sigma))
    @staticmethod
    def pdf(mu,sigma,x):
        return 1/(x*sigma*math.sqrt(2*math.pi))*math.exp(-(math.log(x)-mu)**2/(2*sigma**2))
    @staticmethod
    def cdf(mu,sigma,x):
        return 1/2+1/2*math.erf((math.log(x)-mu)/math.sqrt(2)*sigma)
    @staticmethod
    def kurtosis(mu,sigma):
        return math.exp(4*sigma**2)+2*math.exp(3*sigma**2)+3*math.exp(2*sigma**2)-6
    @staticmethod
    def mean(mu,sigma):
        return math.exp(mu+sigma**2/2)
    @staticmethod
    def median(mu,sigma):
        return math.exp(mu)
    @staticmethod
    def mode(mu,sigma):
        return math.exp(mu-sigma**2)
    @staticmethod
    def variance(mu,sigma):
        return (math.exp(sigma**2)-1)*math.exp(2*mu+sigma**2)
    @staticmethod
    def stddev(mu,sigma):
        return math.sqrt(math.exp(sigma**2)-1)*math.exp(2*mu+sigma**2)
    @staticmethod
    def entropy(mu,sigma):
        return math.log(sigma*math.exp(mu+1/2)*math.sqrt(2*math.pi))
    @staticmethod
    def skewness(mu,sigma):
        return (math.exp(sigma**2)+2)*math.sqrt(math.exp(sigma**2)-1)
    @staticmethod
    def mle(x):
        params={'mu':np.sum(np.log(np.array(x)))/len(x)}
        params['sigma2']=np.sum(np.power(np.log(np.array(x))-params['mu'],2))
        return params
