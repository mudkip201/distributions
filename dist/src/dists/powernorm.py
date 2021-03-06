'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class powernorm(Distribution): #power-normal
    @staticmethod
    def random(c):
        return st.powernorm.rvs(c)