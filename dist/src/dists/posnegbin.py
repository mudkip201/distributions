'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.negbin.negbin as negbin

class posnegbin(Distribution): #positive negative binomial (returns only >0)
    @staticmethod
    def random(k,p):
        n=negbin.random(k,p)
        while(n==0):
            n=negbin.random(k,p)
        return n
