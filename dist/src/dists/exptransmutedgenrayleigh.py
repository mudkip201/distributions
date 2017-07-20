'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class exptransmutedgenrayleigh(Distribution): #exponential transmuted generalized rayleigh
    @staticmethod
    def random(a,b,d,l):
        i=1+l-math.sqrt((1+l)**2-4*l*math.pow(ds.rg0(),1/d))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a)))/b
    @staticmethod
    def pdf(a,b,d,l,x):
        2*a*d*b**2*x*math.exp(-(b*x)**2)*math.pow(1-math.exp(-(b*x)**2),a*d-1)*(1+l-2*l*math.pow(1-math.exp(-(b*x)**2),a))*math.pow(1+l-l*math.pow(1-math.exp(-(b*x)**2),a),d-1)
    @staticmethod
    def cdf(a,b,d,l,x):
        return math.pow(1-math.exp(-(b*x)**2),a*d)*math.pow(1+l-l*math.pow(1-math.exp(-(b*x)**2),a),d)
    @staticmethod
    def median(a,b,d,l):
        i=1+l-math.sqrt((1+l)**2-4*l*math.pow(1/2,1/d))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a)))/b
    @staticmethod
    def ppf(a,b,d,l,q):
        i=1+l-math.sqrt((1+l)**2-4*l*math.pow(q,1/d))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a)))/b
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=exptransmutedgenrayleigh.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'d':ret[2],'l':ret[3]}
