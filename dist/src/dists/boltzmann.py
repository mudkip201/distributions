'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class boltzmann(Distribution):
    @staticmethod
    def random(N,l):
        return math.ceil(-math.log(1-ds.rg0()*(1-math.exp(-l*N)))/l-1)
    @staticmethod
    def pdf(N,l,x):
        return (1-math.exp(-l))/(1-math.exp(-l*N))*math.exp(-l*x)
    @staticmethod
    def cdf(N,l,x):
        if(x>=N-1):
            return 1
        return (1-math.exp(-l*(math.floor(x)+1)))/(1-math.exp(-l*N))
    @staticmethod
    def median(N,l):
        return math.ceil(-math.log(1-(1-math.exp(-l*N))/2)/l-1)
    @staticmethod
    def ppf(N,l,q):
        return math.ceil(-math.log(1-q*(1-math.exp(-l*N)))/l-1)
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=boltzmann.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(-1,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(-10,-0.001),(0.001,10)]).x.tolist()
        return {'N':ret[0],'l':ret[1]}