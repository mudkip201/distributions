'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class pearson3(Distribution):
    @staticmethod
    def random(skew):
        return st.pearson3.rvs(skew)