'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''
import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.stats as st
import scipy.optimize as op

'''
Alpha distribution
'''
class alpha(Distribution):
    @staticmethod
    def random(aa):
        n=ds.rg0()
        return 1/(aa-st.norm.ppf(n*st.norm.cdf(aa)))
    @staticmethod
    def pdf(aa,x):
        return 1/(x**2*st.norm.cdf(aa)*math.sqrt(2*math.pi))*math.exp(-1/2*(aa-1/x)**2)
    @staticmethod
    def cdf(aa,x):
        return st.norm.cdf(aa-1/x)/st.norm.cdf(aa)
    @staticmethod
    def median(aa):
        return 1/(aa-st.norm.ppf(st.norm.cdf(aa)/2))
    @staticmethod
    def ppf(aa,q):
        return 1/(aa-st.norm.ppf(q*st.norm.cdf(aa)))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=alpha.pdf(args_[0],i)
            return -tomin
        ret=op.minimize(mlefunc,(1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0]}