'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class modburr3kumaraswamy(Distribution):
    @staticmethod
    def random(a,aa,b,bb,g):
        u=ds.rg0()
        return math.pow(1-math.pow(math.pow((math.pow(u,-g/aa)-1)/g,-1/bb)+1,-1/b),1/a)
    @staticmethod
    def pdf(a,aa,b,bb,g,x):
        r=math.pow(1-math.pow(x,a),-b)-1
        return aa*bb*a*b*math.pow(x,a-1)*math.pow(1-math.pow(x,a),-b-1)*math.pow(r,-bb-1)*math.pow(1+g*math.pow(r,-bb),-aa/g-1)
    @staticmethod
    def cdf(a,aa,b,bb,g,x):
        return math.pow(1+g*math.pow(math.pow(1-math.pow(x,a),-b)-1,-bb),-aa/g)
    @staticmethod
    def median(a,aa,b,bb,g):
        return math.pow(1-math.pow(math.pow((math.pow(1/2,-g/aa)-1)/g,-1/bb)+1,-1/b),1/a)
    @staticmethod
    def ppf(a,aa,b,bb,g,q):
        return math.pow(1-math.pow(math.pow((math.pow(q,-g/aa)-1)/g,-1/bb)+1,-1/b),1/a)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=modburr3kumaraswamy.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'aa':ret[1],'b':ret[2],'bb':ret[3],'g':ret[4]}