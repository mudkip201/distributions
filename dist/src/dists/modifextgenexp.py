'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class modifextgenexp(Distribution): #modified extended generalized exponential
    @staticmethod
    def random(a,b,e,l):
        p=ds.rg0()
        return -1/l*math.log(1/e*(b-math.pow(p*math.pow(b,a/e)-(p-1)*math.pow(b-e,a/e)),e/a))
    @staticmethod
    def pdf(a,b,e,l,x):
        return a*l*math.pow(b-e*math.exp(-l*x),a/e-1)*math.exp(-l*x)/(math.pow(b,a/e)-math.pow(b-e,a/e))
    @staticmethod
    def cdf(a,b,e,l,x):
        return (math.pow(b-e*math.exp(-l*x),a/e)-math.pow(b-e,a/e))/(math.pow(b,a/e)-math.pow(b-e,a/e))
    @staticmethod
    def median(a,b,e,l):
        return -1/l*math.log((1/e)*(b-math.pow((math.pow(b,a/e)+math.pow(b+e,a/e))/2,e/a)))
    @staticmethod
    def mode(a,b,e,l):
        if(a>e):
            return -1/l*math.log(b/a)
    @staticmethod
    def ppf(a,b,e,l,q):
        return -1/l*math.log(1/e*(b-math.pow(q*math.pow(b,a/e)-(q-1)*math.pow(b-e,a/e)),e/a))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=modifextgenexp.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'e':ret[2],'l':ret[3]}