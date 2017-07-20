'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.normal.normal as normal

class chi2(Distribution): #Chi-squared
    @staticmethod
    def random(k):
        avg_=0
        for _ in range(k):
            avg_+=math.pow(normal.random(0,1),2)
        return avg_
    @staticmethod
    def pdf(k,x):
        return 1/(math.pow(2,k/2)*math.gamma(k/2))*math.pow(x,k/2-1)*math.exp(-x/2)
    @staticmethod
    def kurtosis(k):
        return 12/k
    @staticmethod
    def mean(k):
        return k
    @staticmethod
    def median(k):
        return k*math.pow(1-2/(9*k),3)
    @staticmethod
    def mode(k):
        return max(k-2,0)
    @staticmethod
    def variance(k):
        return 2*k
    @staticmethod
    def stddev(k):
        return math.sqrt(2*k)
    @staticmethod
    def skewness(k):
        return math.sqrt(8/k)
    @staticmethod
    def mle(x): #not in eclipse
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=chi2.pdf(args_[0],i)
            return -tomin
        '''
        ret=op.differential_evolution(mlefunc,[(0.01,50)]).x.tolist()
        return {'k':ret[0]}
'''