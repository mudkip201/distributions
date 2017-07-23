'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class reflectedtoppleone(Distribution):
    @staticmethod
    def random(aa,bb):
        n=ds.rg0()
        if aa==1:
            return 1-math.pow(1-n,1/bb)
        if(0<aa and aa<1):
            return 1-(aa+math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1-n,1/bb)))/(2*(aa-1))
        return 1-(aa-math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1-n,1/bb)))/(2*(aa-1))
    @staticmethod
    def median(aa,bb):
        if aa==1:
            return 1-math.pow(1/2,1/bb)
        if(0<aa and aa<1):
            return 1-(aa+math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1/2,1/bb)))/(2*(aa-1))
        return 1-(aa-math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1/2,1/bb)))/(2*(aa-1))
    @staticmethod
    def ppf(aa,bb,q):
        if aa==1:
            return 1-math.pow(1-q,1/bb)
        if(0<aa and aa<1):
            return 1-(aa+math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1-q,1/bb)))/(2*(aa-1))
        return 1-(aa-math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1-q,1/bb)))/(2*(aa-1))
