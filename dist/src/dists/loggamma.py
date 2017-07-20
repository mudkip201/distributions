'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class loggamma(Distribution):
    @staticmethod
    def random(c):
        return st.loggamma.rvs(c)
