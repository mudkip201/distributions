'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class kumaraswamygenexppareto(Distribution): #kumaraswamy generalized exponentiated pareto
    @staticmethod
    def random(a,b,l,t):
        return math.pow(1-math.pow(1-math.pow(1-ds.rg0(),1/b),1/(a*t)),-1/l)-1
    @staticmethod
    def pdf(a,b,l,t,x):
        return a*b*t*l*math.pow(1-math.pow(1+x,-l),t*a-1)*math.pow(1+x,-l-1)*math.pow(1-math.pow(1-math.pow(1+x,-l),t*a),b-1)
    @staticmethod
    def cdf(a,b,l,t,x):
        return 1-math.pow(1-math.pow(1-math.pow(1+x,-l),t*a),b)
    @staticmethod
    def median(a,b,l,t):
        return math.pow(1-math.pow(1-math.pow(1/2,1/b),1/(a*t)),-1/l)-1
    @staticmethod
    def ppf(a,b,l,t,q):
        return math.pow(1-math.pow(1-math.pow(1-q,1/b),1/(a*t)),-1/l)-1
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamygenexppareto.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2],'t':ret[3]}