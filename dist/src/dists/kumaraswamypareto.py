'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class kumaraswamypareto(Distribution):
    @staticmethod
    def random(a,b,bb,k):
        return bb/math.pow(1-math.pow(1-math.pow(1-ds.rg0(),1/b),1/a),1/k)
    @staticmethod
    def pdf(a,b,bb,k,x):
        return a*b*k*math.pow(bb,k)/math.pow(x,k+1)*math.pow(1-math.pow(bb/x,k),a-1)*math.pow(1-math.pow(1-math.pow(bb/x,k),a),b-1)
    @staticmethod
    def cdf(a,b,bb,k,x):
        return 1-math.pow(1-math.pow(1-math.pow(bb/x,k),a),b)
    @staticmethod
    def median(a,b,bb,k):
        return bb/math.pow(1-math.pow(1-math.pow(1/2,1/b),1/a),1/k)
    @staticmethod
    def ppf(a,b,bb,k,q):
        return bb/math.pow(1-math.pow(1-math.pow(1-q,1/b),1/a),1/k)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamypareto.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'bb':ret[2],'k':ret[3]}
