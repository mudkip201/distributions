'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class exponentiatedexpbin(Distribution): #exponentiated exponential binomial
    @staticmethod
    def random(a,l,n,t):
        return -1/l*math.log(1-math.pow(1/t*(1-math.pow(1-(1-math.pow(1-t,n))*ds.rg0(),1/n)),1/a))
    @staticmethod
    def pdf(a,l,n,t,x):
        return l*a*n*t*math.exp(-l*x)*math.pow(1-math.exp(-l*x),a-1)*math.pow(1-t*math.pow(1-math.exp(-l*x),a),n-1)/(1-math.pow(1-t,n))
    @staticmethod
    def cdf(a,l,n,t,x):
        return (1-math.pow(1-t*math.pow(1-math.exp(-l*x),a),n))/(1-math.pow(1-t,n))
    @staticmethod
    def median(a,l,n,t):
        return -1/l*math.log(1-math.pow(1/t*(1-math.pow(1-(1-math.pow(1-t,n))/2,1/n)),1/a))
    @staticmethod
    def ppf(a,l,n,t,q):
        return -1/l*math.log(1-math.pow(1/t*(1-math.pow(1-(1-math.pow(1-t,n))*q,1/n)),1/a))
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=exponentiatedexpbin.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,0.5,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'l':ret[1],'n':ret[2],'t':ret[3]}
