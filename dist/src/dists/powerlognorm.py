'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.stats as st


class powerlognorm(Distribution): #power-log normal
    @staticmethod
    def random(c,s):
        return st.powerlognorm.rvs(c,s)