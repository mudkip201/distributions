'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.negbin.negbin as negbin
import dists.beta.beta as beta


class betapascal(Distribution):
    @staticmethod
    def random(k,m,n):
        return negbin.random(k,beta.random(m,n))