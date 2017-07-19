'''
Created on Jul 16, 2017

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

class transmutedexponentiatedlomax(Distribution): #transmuted exponentiated lomax
    @staticmethod
    def random(a,g,l,t):
        u=ds.rg0()
        i=(1+l-math.sqrt((1+l)**2-4*l*u))/(2*l)
        return (math.pow(1-math.pow(i,1/a),-1/t)-1)/g
    @staticmethod
    def pdf(a,g,l,t,x):
        return a*t*g*math.pow(1-math.pow(1+g*x,-t),a-1)/(math.pow(1+g*x,t+1))*(1+l-2*l*math.pow(1-math.pow(1+g*x,-t),a))
    @staticmethod
    def cdf(a,g,l,t,x):
        return math.pow(1-math.pow(1+g*x,-t),a)*(1+l-l*math.pow(1-math.pow(1+g*x,-t),a))
    @staticmethod
    def mean(a,g,l,t):
        return a*(1+l)/g*(sp.beta(1-1/t,a)-1/a)-2*a*l/g*(sp.beta(1-1/t,2*a)-1/(2*a))
    @staticmethod
    def median(a,g,l,t):
        i=(1+l-math.sqrt((1+l)**2-2*l))/(2*l)
        return (math.pow(1-math.pow(i,1/a),-1/t)-1)/g
    @staticmethod
    def ppf(a,g,l,t,q):
        i=(1+l-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return (math.pow(1-math.pow(i,1/a),-1/t)-1)/g
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=transmutedexponentiatedlomax.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'g':ret[1],'l':ret[2],'t':ret[3]}