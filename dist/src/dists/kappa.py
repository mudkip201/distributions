'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class kappa(Distribution):
    @staticmethod
    def random(h,k,xi,aa):
        n=ds.rg0()
        return xi+(aa/k)*(1-math.pow(1-math.pow(n,h),k))
    @staticmethod
    def median(h,k,xi,aa):
        return xi+(aa/k)*(1-math.pow(1-math.pow(1/2,h),k))
    @staticmethod
    def ppf(h,k,xi,aa,q):
        return xi+(aa/k)*(1-math.pow(1-math.pow(q,h),k))
