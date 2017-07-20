'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class expgeninvweibull(Distribution): #exponentiated generalized inverse weibull
    @staticmethod
    def random(a,b,l,t):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(ds.rg0(),1/b),1/a)),-1/t)
    @staticmethod
    def pdf(a,b,l,t,x):
        return a*b*t*math.pow(l,t)*math.pow(x,-t-1)*math.exp(-math.pow(l/x,t))*math.pow(1-math.exp(-math.pow(l/x,t)),a-1)*math.pow(1-math.pow(1-math.exp(math.pow(l/x,t)),a),b-1)
    @staticmethod
    def cdf(a,b,l,t,x):
        return math.pow(1-math.pow(1-math.exp(math.pow(l/x,t)),a),b)
    @staticmethod
    def median(a,b,l,t):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(1/2,1/b),1/a)),-1/t)
    @staticmethod
    def ppf(a,b,l,t,q):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(q,1/b),1/a)),-1/t)
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expgeninvweibull.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2],'t':ret[3]}