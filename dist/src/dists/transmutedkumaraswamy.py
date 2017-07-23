'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class transmutedkumaraswamy(Distribution): #transmuted kumaraswamy
    @staticmethod
    def random(a,l,t):
        q=ds.rg0()
        return math.pow(1-math.pow(1-((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l),1/t),1/a)
    @staticmethod
    def pdf(a,l,t,x):
        return a*t*math.pow(x,a-1)*math.pow(1-math.pow(x,a),t-1)*(1-l+2*l*math.pow(1-math.pow(x,a),t))
    @staticmethod
    def cdf(a,l,t,x):
        return (1-math.pow(1-math.pow(x,a),t))*(1+l*math.pow(1-math.pow(x,a),t))
    @staticmethod
    def median(a,l,t):
        return math.pow(1-math.pow(1-((1+l)-math.sqrt(1+l**2))/(2*l),1/t),1/a)
    @staticmethod
    def ppf(a,l,t,q):
        return math.pow(1-math.pow(1-((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l),1/t),1/a)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=transmutedkumaraswamy.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'l':ret[1],'t':ret[2]}