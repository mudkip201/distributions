'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class gengamma(Distribution): #generalized gamma
    @staticmethod
    def random(a,c):
        return st.gengamma.rvs(a,c)