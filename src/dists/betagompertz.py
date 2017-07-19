'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import numpy as np
from numpy import random as r
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as op
import dists.beta.beta as beta

class betagompertz(Distribution): 
    @staticmethod
    def random(g,t,a,b):
        return (1/g)*math.log(1-g/t*math.log(1-beta.random(a,b)))
    @staticmethod
    def pdf(g,t,a,b,x):
            return t*math.exp(g*x)*math.exp(-b*t/g*(math.exp(g*x)-1))/sp.beta(a,b)*math.pow(1-math.exp(-t/g*(math.exp(g*x)-1)),a-1)
    @staticmethod
    def mle(x): #better, but still not good
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=betagompertz.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(0.25,0.25,0.25,0.25),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.001,10),(0.001,10),(0.001,10),(0.001,10)]).x.tolist()
        return {'g':ret[0],'t':ret[1],'a':ret[2],'l':ret[3]}