'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.negbin.negbin as negbin

class bhattacharyanegbin(Distribution): #bhattacharya negative binomial
    @staticmethod
    def random(a,b,k):
        return negbin.random(k+negbin.random(a,b/(b+1)),(b+1)/(b+2))