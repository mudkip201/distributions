'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class hypergeo(Distribution): #hypergeometric
    @staticmethod
    def random(n,N,K):
        num_success=0
        for i in range(n):
            d=r.random()
            if(d<(K-num_success)/(N-i)):
                num_success+=1
        return num_success
    @staticmethod
    def kurtosis(n,N,K):
        return 1/(n*K*(N-K)*(N-n)*(N-2)*(N-3))*((N-1)*N**2*(N*(N+1)-6*K*(N-K)-6*n*(N-n))+6*n*K*(N-K)*(N-n)*(5*N-6))
    @staticmethod
    def mean(n,N,K):
        return n*K/N
    @staticmethod
    def mode(n,N,K):
        return math.floor((n+1)*(K+1)/(N+2))
    @staticmethod
    def variance(n,N,K):
        return n*K/N*(N-K)/N*(N-n)/(N-1)
    @staticmethod
    def stddev(n,N,K):
        return math.sqrt(n*K/N*(N-K)/N*(N-n)/(N-1))
    @staticmethod
    def skewness(n,N,K):
        return (N-2*K)*math.sqrt(N-1)*(N-2*n)/(math.sqrt(n*K*(N-K)*(N-n))*(N-2))