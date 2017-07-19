'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import numpy as np
from numpy import random as r
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as op

class neghypergeo(Distribution): #negative hypergeometric
    @staticmethod
    def random(R,N,K):
        num_success=0
        num_failure=0
        total_=0
        while(num_failure<R):
            if(r.random()<=(K-num_success)/(N-total_)):
                num_success+=1
            else:
                num_failure+=1
                total_+=1
        return num_success
    @staticmethod
    def mean(R,N,K):
        return R*K/(N-K+1)
    @staticmethod
    def variance(R,N,K):
        return R*(N+1)*K/((N-K+1)*(N-K+2))*(1-R/(N-K+1))
    @staticmethod
    def stddev(R,N,K):
        return math.sqrt(R*(N+1)*K/((N-K+1)*(N-K+2))*(1-R/(N-K+1)))