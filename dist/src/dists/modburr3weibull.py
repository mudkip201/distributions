'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class modburr3weibull(Distribution):
    @staticmethod
    def random(a,b,g,l,n):
        return math.pow(math.log(math.pow(math.pow(ds.rg0(),-g/a)-1,-1/b)/g+1)/l,1/n)
    @staticmethod
    def pdf(a,b,g,l,n,x):
        return a*b*math.pow(math.exp(l*math.pow(x,n)-1),-b-1)*l*n*math.pow(x,n-1)*math.exp(l*math.pow(x,n))*math.pow(1+g*math.pow(math.exp(l*math.pow(x,n))-1,-b),-a/g-1)
    @staticmethod
    def median(a,b,g,l,n):
        return math.pow(math.log(math.pow(math.pow(1/2,-g/a)-1,-1/b)/g+1)/l,1/n)
    @staticmethod
    def ppf(a,b,g,l,n,q):
        return math.pow(math.log(math.pow(math.pow(q,-g/a)-1,-1/b)/g+1)/l,1/n)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=modburr3weibull.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'g':ret[2],'l':ret[3],'n':ret[4]}