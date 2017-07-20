'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class asymlaplace(Distribution): #asymmetric laplace
    @staticmethod
    def pdf(m,lmbda,kppa,x):
        if(kppa<=0 or lmbda<=0):
            raise ValueError("kppa and lambda must be greater than 0")
        if(x>=m):
            return lmbda/(kppa+1/kppa)*math.exp(-lmbda*kppa*(x-m))
        return lmbda/(kppa+1/kppa)*math.exp((lmbda/kppa)*(x-m))
    @staticmethod
    def cdf(m,lmbda,kppa,x):
        if(kppa<=0 or lmbda<=0):
            raise ValueError("kppa and lambda must be greater than 0")
        if(x>m):
            return 1-1/(1+kppa**2)*math.exp(-lmbda*kppa*(x-m))
        return kppa**2/(1+kppa**2)*math.exp(lmbda/kppa*(x-m))
    @staticmethod
    def random(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise ValueError("kppa and lambda must be greater than 0")
        n=ds.rg0()
        n=n*(1/kppa+kppa)-kppa
        s=abs(n)/n
        return m-1/(lmbda*s*math.pow(kppa,s))*math.log(n*s*math.pow(kppa,s))
    @staticmethod
    def mean(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise ValueError("kppa and lambda must be greater than 0")
        return m+(1-kppa**2)/(lmbda*kppa)
    @staticmethod
    def median(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise ValueError("kppa and lambda must be greater than 0")
        return m+(kppa/lmbda)*math.log((1+kppa**2)/(2*kppa**2))
    @staticmethod
    def variance(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise ValueError("kppa and lambda must be greater than 0")
        return (1+math.pow(kppa,4))/(lmbda**2*kppa**2)
    @staticmethod
    def stddev(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise ValueError("kppa and lambda must be greater than 0")
        return math.sqrt((1+math.pow(kppa,4))/(lmbda**2*kppa**2))
    @staticmethod
    def skewness(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise ValueError("kppa and lambda must be greater than 0")
        return 2*(1-math.pow(kppa,6))/(math.pow(math.pow(kppa,4)+1,3/2))
    @staticmethod
    def entropy(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise ValueError("kppa and lambda must be greater than 0")
        return math.log(math.exp(1)*(1+kppa**2)/(kppa*lmbda))
    #@staticmethod
    #def ppf(m,lmbda,kppa,q): #not working
        #if(kppa<=0 or lmbda<=0):
            #raise ValueError("kppa and lambda must be greater than 0")
        #if(q<=(math.pow(kppa,2))/(1+math.pow(kppa,2))):
            #return kppa/math.sqrt(2)*math.log((1+math.pow(kppa,2))/math.pow(k,2)*q)
        #return 1/(kppa*math.sqrt(2))*math.log((1+math.pow(kppa,2))*(1-q))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=asymlaplace.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(0,1,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(-10,10),(0.001,10),(0.001,10)]).x.tolist()
        return {'m':ret[0],'lambda':ret[1],'kappa':ret[2]}