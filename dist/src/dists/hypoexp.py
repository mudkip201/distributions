'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.exponential.exponential as exponential

class hypoexp(Distribution): #hypoexponential
    @staticmethod
    def random(lmbdas):
        sum_=0
        for i in range(lmbdas):
            sum_+=exponential.random(lmbdas[i])
        return sum_
