'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class invweibull(Distribution): #inverse weibull
    @staticmethod
    def random(c):
        return st.invweibull.rvs(c)
