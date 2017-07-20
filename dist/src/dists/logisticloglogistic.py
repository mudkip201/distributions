'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import scipy.optimize as op
import dists.logisticburr12.logisticburr12 as logisticburr12

class logisticloglogistic(Distribution):
    @staticmethod
    def random(c,l,s):
        return logisticburr12.random(c,1,l,s)
    @staticmethod
    def pdf(c,l,s,x):
        return logisticburr12.pdf(c,1,l,s,x)
    @staticmethod
    def cdf(c,l,s,x):
        return logisticburr12.cdf(c,1,l,s,x)
    @staticmethod
    def median(c,l,s):
        return logisticburr12.median(c,1,l,s)
    @staticmethod
    def ppf(c,l,s,q):
        return logisticburr12.ppf(c,1,l,s,q)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logisticloglogistic.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'c':ret[0],'l':ret[1],'s':ret[2]}