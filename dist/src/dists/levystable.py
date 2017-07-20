'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class levystable(Distribution):
    @staticmethod
    def random(a,b):
        return st.levy_stable.rvs(a,b)