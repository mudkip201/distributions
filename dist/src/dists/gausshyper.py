'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st

class gausshyper(Distribution):
    @staticmethod
    def random(a,b,c,z):
        return st.gausshyper.rvs(a,b,c,z)
