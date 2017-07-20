'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class levyl(Distribution):
    @staticmethod
    def random():
        return st.levy_l.rvs()
