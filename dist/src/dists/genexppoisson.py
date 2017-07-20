'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class genexppoisson(Distribution): #generalized exponentiated poisson
    @staticmethod
    def random(a,b,l):
        return -math.log(1+math.log(1-math.pow(ds.rg0(),1/a)*(1-math.exp(-l)))/l)/b
    @staticmethod
    def pdf(a,b,l,x):
        return a*l*b/math.pow(1-math.exp(-l),a)*math.pow(1-math.exp(-l+l*math.exp(-b*x)),a-1)*math.exp(-l-b*x+l*math.exp(-b*x))
    @staticmethod
    def cdf(a,b,l,x):
        return math.pow((1-math.exp(-l+l*math.exp(-b*x)))/(1-math.exp(-l)),a)
    @staticmethod
    def median(a,b,l):
        return -math.log(1+math.log(1-math.pow(1/2,1/a)*(1-math.exp(-l)))/l)/b
    @staticmethod
    def ppf(a,b,l,q):
        return -math.log(1+math.log(1-math.pow(q,1/a)*(1-math.exp(-l)))/l)/b
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=genexppoisson.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2]}