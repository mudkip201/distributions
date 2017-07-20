'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class genlogistic5(Distribution): #generalized logistic V
    @staticmethod
    def random(aa,mu,sigma):
        return mu+sigma*(1/aa)*(1-math.pow(1/ds.rg0()-1,aa))
    @staticmethod
    def pdf(aa,mu,sigma,x):
        z=(x-mu)/sigma
        return math.pow(1-aa*z,(1/aa)-1)/(1+math.pow(1-aa*z,1/aa))**2
    @staticmethod
    def cdf(aa,mu,sigma,x):
        z=(x-mu)/sigma
        return 1/(1+math.pow(1-aa*z,1/aa))
    @staticmethod
    def median(aa,mu,sigma):
        return mu+sigma*(1/aa)*(1-math.pow(1,aa))
    @staticmethod
    def ppf(aa,mu,sigma,q):
        return mu+sigma*(1/aa)*(1-math.pow(1/q-1,aa))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=genlogistic5.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0],'mu':ret[1],'sigma':ret[2]}