'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.pareto4.pareto4 as pareto4

class pareto3(Distribution):
    @staticmethod
    def cdf(mu,sigma,gmma,x):
        return pareto4.cdf(mu,sigma,gmma,1,x)
