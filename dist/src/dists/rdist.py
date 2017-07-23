'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class rdist(Distribution):
    @staticmethod
    def random(c):
        return st.rdist.rvs(c)