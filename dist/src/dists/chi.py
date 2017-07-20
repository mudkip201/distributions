'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.optimize as op
import dists.chi2.chi2 as chi2

class chi(Distribution):
    @staticmethod
    def random(k):
        return math.sqrt(chi2.random(k))
    @staticmethod
    def pdf(k,x):
        return math.pow(2,1-k/2)*math.pow(x,k-1)*math.exp(-x**2/2)/math.gamma(k/2)
    @staticmethod
    def mean(k):
        return math.sqrt(2)*math.gamma((k+1)/2)/math.gamma(k/2)
    @staticmethod
    def mode(k):
        if(k>=1):
            return math.sqrt(k-1)
    @staticmethod
    def variance(k):
        return k-(math.sqrt(2)*math.gamma((k+1)/2)/math.gamma(k/2))**2
    @staticmethod
    def stddev(k):
        return math.sqrt(k-(math.sqrt(2)*math.gamma((k+1)/2)/math.gamma(k/2))**2)
    @staticmethod
    def skewness(k):
        return (k-(math.sqrt(2)*math.gamma((k+1)/2)/math.gamma(k/2))**2)/chi.stddev(k)**3*(1-2*chi.variance(k))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=chi.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,50)]).x.tolist()
        return {'k':ret[0]}