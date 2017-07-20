'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class kolmogovsmirnovtwolarge(Distribution):
    @staticmethod
    def random(n):
        return st.kstwobign.rvs(n)