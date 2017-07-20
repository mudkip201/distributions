'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class expnorm(Distribution): #exponential normal
    @staticmethod
    def random(k):
        return st.exponnorm.rvs(k)
