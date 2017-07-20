'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.optimize as op
import dists.normal.normal as normal

class johnsonsb(Distribution):
    @staticmethod
    def random(delta,gmma,xi,lmbda):
        v=(normal.random(0,1)-gmma)/delta
        return xi+lmbda*(1/(1+math.exp(-v)))
    @staticmethod
    def pdf(delta,gmma,xi,lmbda,x):
        return delta*lmbda*math.exp(-1/2*(gmma+delta*math.log((x-xi)/(-x+xi+lmbda)))**2)/(math.sqrt(2*math.pi)*(x-xi)*(-x+xi+lmbda))
    @staticmethod
    def cdf(delta,gmma,xi,lmbda,x):
        if(xi<x and x<xi+lmbda/2):
            return 1/2*math.erfc(-(gmma+delta*math.log((x-xi)/(-x+xi+lmbda)))/math.sqrt(2))
        if(xi+lmbda/2<=x and x<xi+lmbda):
            return 1/2*(1+math.erf((gmma+delta*math.log((x-xi)/(-x+xi+lmbda)))/math.sqrt(2)))
    @staticmethod
    def median(delta,gmma,xi,lmbda):
        return xi+lmbda/(1+math.exp(gmma/delta))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=johnsonsb.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'delta':ret[0],'gamma':ret[1],'xi':ret[2],'lmbda':ret[3]}