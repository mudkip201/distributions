'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.gamma.gamma as gamma

class scaledchi2(Distribution): #scaled chi-squared
    @staticmethod
    def random(s,k):
        return gamma.random(2*(s**2),k/2.0)