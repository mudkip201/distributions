'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class kolmogorovsmirnovone(Distribution):
    @staticmethod
    def random(n):
        return st.ksone.rvs(n)
