'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class kumaraswamygenpowerweibull(Distribution):#kumaraswamy generalized power-weibull
    @staticmethod
    def random(a,b,t,aa,l):
        return l*math.pow(math.pow(1-math.log(1-math.pow(1-math.pow(1-ds.rg0(),1/b),1/a)),1/t)-1,1/aa)
    @staticmethod
    def pdf(a,b,t,aa,l,x):
        p5=math.pow(1-math.pow(1-math.exp(1-math.pow(1+math.pow(x/l,aa),t)),a),b-1)
        p4=math.pow(1-math.exp(1-math.pow(1+math.pow(x/l,aa),t)),a-1)
        p3=math.exp(1-math.pow(1+math.pow(x/l,aa),t))
        p2=math.pow(1+math.pow(x/l,aa),t-1)
        p1=a*b*aa*math.pow(x,aa-1)/math.pow(l,aa)
        return p1*p2*p3*p4*p5
    @staticmethod
    def cdf(a,b,t,aa,l,x):
        return 1-math.pow(1-math.pow(1-math.exp(1-math.pow(1+math.pow(x/l,aa),t)),a),b)
    @staticmethod
    def median(a,b,t,aa,l):
        return l*math.pow(math.pow(1-math.log(1-math.pow(1-math.pow(1/2,1/b),1/a)),1/t)-1,1/aa)
    @staticmethod
    def ppf(a,b,t,aa,l,q):
        return l*math.pow(math.pow(1-math.log(1-math.pow(1-math.pow(1-q,1/b),1/a)),1/t)-1,1/aa)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamygenpowerweibull.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'t':ret[2],'aa':ret[3],'l':ret[4]}