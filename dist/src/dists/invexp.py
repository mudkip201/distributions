'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.invgamma.invgamma as invgamma

class invexp(Distribution): #inverse exponential
    @staticmethod
    def random(theta):
        return invgamma.random(theta,1)