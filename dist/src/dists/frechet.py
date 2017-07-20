'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class frechet(Distribution):
    @staticmethod
    def pdf(aa,m,s,x):
        if(aa<=0 or s<=0):
            raise ValueError("aa and s must be positive")
        if(x<=m):
            raise ValueError("x must be bigger than m")
        return aa/s*math.pow((x-m)/s,-1-aa)*math.exp(-math.pow((x-m)/s,-aa))
    @staticmethod
    def cdf(aa,m,s,x):
        if(aa<=0 or s<=0):
            raise ValueError("aa and s must be positive")
        if(x<=m):
            raise ValueError("x must be bigger than m")
        return math.exp(-math.pow((x-m)/s,-aa))
    @staticmethod
    def random(aa,m,s):
        if(aa<=0 or s<=0):
            raise ValueError("aa and s must be positive")
        n=ds.rg0()
        return -math.pow(math.log(n),-1/aa)*(m*math.pow(-math.log(n),1/aa)+s)
    @staticmethod
    def median(aa,m,s):
        if(aa<=0 or s<=0):
            raise ValueError("aa and s must be positive")
        return m+s/math.pow(math.log(2),1/aa)
    @staticmethod
    def mode(aa,m,s):
        if(aa<=0 or s<=0):
            raise ValueError("aa and s must be positive")
        return m+s*math.pow(aa/(1+aa),1/aa)
    @staticmethod
    def entropy(aa,m,s):
        if(aa<=0 or s<=0):
            raise ValueError("aa and s must be positive")
        return 1+ds.euler_gamma/aa+ds.euler_gamma+math.log(s/aa)
    @staticmethod
    def ppf(aa,m,s,q):
        if(aa<=0 or s<=0):
            raise ValueError("aa and s must be positive")
        return -math.pow(math.log(q),-1/aa)*(m*math.pow(-math.log(q),1/aa)+s)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=frechet.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0],'m':ret[1],'s':ret[2]}