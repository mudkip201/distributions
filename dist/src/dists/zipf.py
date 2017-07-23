'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class zipf(Distribution):
    @staticmethod
    def random(a):
        return st.zipf.rvs(a)
