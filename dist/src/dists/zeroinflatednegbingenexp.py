'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import dists.negbin.negbin as negbin

class zeroinflatednegbingenexp(Distribution): #zero-inflated negative binomial-generalized exponential
    @staticmethod
    def random(a,b,f,r):
        u=ds.rg0()
        l=-1/b*math.log(1-math.pow(u,1/a))
        y=negbin.random(r,math.exp(-l))
        uu=ds.rg0()
        if(uu>f):
            return y
        return 0
