'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.chi2.chi2 as chi2

class upsilon(Distribution):
    @staticmethod
    def random(ts, nus):
        y=0
        for i in range(len(ts)):
            y+=ts[i]*math.sqrt(chi2.random(nus[i]))
        return y