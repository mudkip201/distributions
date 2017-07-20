'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class expgenfrechet(Distribution): #exponentiated generalized frechet
    @staticmethod
    def random(a,b,l,s):
        return s*l/math.log(1-math.pow(1-math.pow(ds.rg0(),1/b),1/a))
    @staticmethod
    def pdf(a,b,l,s,x):
        return a*b*l*math.pow(s,l)*math.pow(x,-l-1)*math.exp(-math.pow(s/x,l))*math.pow(1-math.exp(-math.pow(s/x,l)),a-1)*math.pow(1-math.pow(1-math.exp(-math.pow(s/x,l)),a),b-1)
    @staticmethod
    def cdf(a,b,l,s,x):
        return math.pow(1-math.pow(1-math.exp(-math.pow(s/x,l)),a),b)
    @staticmethod
    def median(a,b,l,s):
        return s*l/math.log(1-math.pow(1-math.pow(1/2,1/b),1/a))
    @staticmethod
    def ppf(a,b,l,s,q):
        return s*l/math.log(1-math.pow(1-math.pow(q,1/b),1/a))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expgenfrechet.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,5),(0.01,5),(0.01,5),(0.01,5)])
        return {'a':ret[0],'b':ret[1],'l':ret[2],'s':ret[3]}