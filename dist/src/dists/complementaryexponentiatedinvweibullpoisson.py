'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class complementaryexponentiatedinvweibullpoisson(Distribution):
    @staticmethod
    def random(b,l,t):
        return math.pow(-math.log(math.log(ds.rg0()*(math.exp(l)-1)+1)/l)/t,-1/b)
    @staticmethod
    def pdf(b,l,t,x): #something wrong here?
        return t*b*l/(math.exp(l)-1)*math.pow(x,-b-1)*math.exp(-math.pow(x,-t*b))*math.exp(l*math.exp(-t*math.pow(x,-b)))
    @staticmethod
    def cdf(b,l,t,x):
        return (math.exp(l*math.exp(-t*math.pow(x,-b)))-1)/(math.exp(l)-1)
    @staticmethod
    def median(b,l,t):
        return math.pow(-math.log(math.log((math.exp(l)-1)/2+1)/l)/t,-1/b)
    @staticmethod
    def ppf(b,l,t,q):
        return math.pow(-math.log(math.log(q*(math.exp(l)-1)+1)/l)/t,-1/b)
    @staticmethod
    def mle(x): #not in eclipse
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=complementaryexponentiatedinvweibullpoisson.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        #ret=op.minimize(mlefunc,(2,2,2),method="Nelder-Mead").x.tolist()
        '''
        ret=op.differential_evolution(mlefunc,[(0.01,10),(0.01,10),(0.01,10)]).x.tolist()
        return {'b':ret[0],'l':ret[1],'t':ret[2]}
        '''