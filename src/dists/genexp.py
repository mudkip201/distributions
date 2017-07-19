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

class genexp(Distribution): #generalized exponential
    @staticmethod
    def random(l,a):
        return -math.log(1-math.pow(ds.rg0(),1/a))/l
    @staticmethod
    def median(l,a):
        return -math.log(1-math.pow(1/2,1/a))/l
    @staticmethod
    def ppf(l,a,q):
        return -math.log(1-math.pow(q,1/a))/l
    @staticmethod
    def mle(x):
        #args_ - (a,l)
        args0=[1,1]
        func= lambda n,x,args_: n*args_[0]+n*args_[1]+(args_[0]-1)*np.sum(np.log(1-np.exp(-args_[1]*x)))-args_[1]*sum(x)
        ret=op.basinhopping(func,args0)
        return {'a':ret[0],'l':ret[1]}