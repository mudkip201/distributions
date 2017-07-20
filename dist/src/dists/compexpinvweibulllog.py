'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class compexpinvweibulllog(Distribution):
    @staticmethod
    def random(b,l,t):
        return math.pow(-math.log((1-(1-l)*ds.rg0())/l)/t,-1/b)
    @staticmethod
    def pdf(b,l,t,x):
        return t*b*l/(-math.log(1-l)*(1-l*math.exp(-t*math.pow(x,-b))))*math.pow(x,-b-1)*math.exp(-t*math.pow(x,-b))
    @staticmethod
    def cdf(b,l,t,x):
        return (math.exp(l*math.exp(-t*math.pow(x,-b)))-1)/(math.exp(l)-1)
    @staticmethod
    def median(b,l,t):
        return math.pow(-math.log((1-(1-l)/2)/l)/t,-1/b)
    @staticmethod
    def ppf(b,l,t,q):
        return math.pow(-math.log((1-(1-l)*q)/l)/t,-1/b)
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=compexpinvweibulllog.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,0.5,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,25),(0.01,0.999),(0.01,25)]).x.tolist()
        return {'b':ret[0],'l':ret[1],'t':ret[2]}