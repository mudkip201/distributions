'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.normal.normal as normal

class noncentralchi2(Distribution): #non-central chi-squared
    @staticmethod
    def random(mu,sigma,nu):
        sum_=0
        for _ in range(nu):
            sum_+=normal.random(mu,sigma)
        return sum_
