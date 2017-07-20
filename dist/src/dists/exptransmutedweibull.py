'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class exptransmutedweibull(Distribution): #exponential transmuted weibull
    @staticmethod
    def random(a,b,l,n):
        i=math.sqrt(1/l*(1-math.pow(ds.rg0(),1/n))+((1-l)/l)**2/4)
        return a*math.pow(-math.log(i-((1-l)/l)/2),1/b)
    @staticmethod
    def pdf(a,b,l,n,x):
        return n*b/a*math.pow(x/a,b-1)*math.exp(-math.pow(x/a,b))*(1-l+2*l*math.exp(-math.pow(x/a,b)))*math.pow(1+(l-1)*math.exp(-math.pow(x/a,b))-l(math.exp(-2(math.pow(x/a,b)))),n-1)
    @staticmethod
    def cdf(a,b,l,n,x):
        return math.pow(1+(l-1)*math.exp(-math.pow(x/a,b))-l(math.exp(-2(math.pow(x/a,b)))),n)
    @staticmethod
    def median(a,b,l,n):
        i=math.sqrt(1/l*(1-math.pow(1/2,1/n))+((1-l)/l)**2/4)
        return a*math.pow(-math.log(i-((1-l)/l)/2),1/b)
    @staticmethod
    def ppf(a,b,l,n,q):
        i=math.sqrt(1/l*(1-math.pow(q,1/n))+((1-l)/l)**2/4)
        return a*math.pow(-math.log(i-((1-l)/l)/2),1/b)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=exptransmutedweibull.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2],'n':ret[3]}