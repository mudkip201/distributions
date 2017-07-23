'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import scipy.special as sp
import scipy.optimize as op

class mcdonaldloglogistic(Distribution):
    @staticmethod
    def random(a,b,c,aa,bb):
        u=r.random()
        while(u==1):
            u=r.random()
        q1=aa*math.pow((u/(1-u)),1/bb) #Q_alpha,beta
        q2=a/c*math.pow((u/(1-u)),1/b) #Q_{a/c},b
        return q1*math.pow(q2,1/c)/(1-q1*math.pow(q2,1/c))
    @staticmethod
    def pdf(a,b,c,aa,bb,x):
        return c/sp.beta(a/c,b)*aa/bb*math.pow(x/bb,a*aa-1)*math.pow(1+math.pow(x/bb,aa),-a-1)*math.pow(1-math.pow(1-1/(1+math.pow(x/bb,aa)),c),b-1)
    @staticmethod
    def median(a,b,c,aa,bb):
        q1=aa #Q_alpha,beta
        q2=a/c #Q_{a/c},b
        return q1*math.pow(q2,1/c)/(1-q1*math.pow(q2,1/c))
    @staticmethod
    def ppf(a,b,c,aa,bb,q):
        q1=aa*math.pow((q/(1-q)),1/bb) #Q_alpha,beta
        q2=a/c*math.pow((q/(1-q)),1/b) #Q_{a/c},b
        return q1*math.pow(q2,1/c)/(1-q1*math.pow(q2,1/c))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=mcdonaldloglogistic.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'c':ret[2],'aa':ret[3],'bb':ret[4]}