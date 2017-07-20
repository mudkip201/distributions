'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.stats as st

class behrensfisher(Distribution):
    @staticmethod
    def random(n1,n2,theta):
        return st.t.rvs(n2)*math.cos(theta)-st.t.rvs(n1)*math.sin(theta)