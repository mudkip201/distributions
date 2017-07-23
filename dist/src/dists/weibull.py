'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class weibull(Distribution):
    @staticmethod
    def pdf(k,lmbda,x):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return k/lmbda*math.pow(x/lmbda,k-1)*math.exp(-math.pow(x/lmbda,k))
    @staticmethod
    def cdf(k,lmbda,x):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        if(x<0):
            raise ValueError("x must be non-negative")
        return 1-math.exp(-math.pow(x/lmbda,k))
    @staticmethod
    def random(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        return lmbda*math.pow((-math.log(ds.rg0())),1/k)
    @staticmethod
    def kurtosis(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        g1=math.gamma(1+1/k)
        g2=math.gamma(1+2/k)
        g3=math.gamma(1+3/k)
        g4=math.gamma(1+4/k)
        return (-6*g1**4+12*g1**2*g2-3*g2**2-4*g1*g3+g4)/((g2-g1**2)**2)
    @staticmethod
    def mean(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        return lmbda*math.gamma(1+1/k)
    @staticmethod
    def median(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        return lmbda*math.pow(math.log(2),1/k)
    @staticmethod
    def mode(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        if(k>1):
            return lmbda*math.pow((k-1)/k,1/k)
        return 0
    @staticmethod
    def variance(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        return lmbda**2*(math.gamma(1+2/k)-(math.gamma(1+1/k)**2))
    @staticmethod
    def stddev(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        return lmbda*math.sqrt(math.gamma(1+2/k)-(math.gamma(1+1/k)**2))
    @staticmethod
    def entropy(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        return ds.euler_gamma*(1-1/k)+math.log(lmbda/k)+1
    @staticmethod
    def skewness(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        return (math.gamma(1+3/k)*lmbda**3-3*weibull.mean(k,lmbda)*weibull.variance(k,lmbda)-weibull.mean(k,lmbda)**3)/math.pow(weibull.variance,3/2)
    @staticmethod
    def ppf(k,lmbda,q):
        if(lmbda<=0 or k<=0):
            raise ValueError("k and lambda must be positive")
        return lmbda*math.pow((-math.log(q)),1/k)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=weibull.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'k':ret[0],'lambda':ret[1]}