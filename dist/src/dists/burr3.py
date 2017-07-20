'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class burr3(Distribution):
    @staticmethod
    def random(r,k,l,s):
        return s*math.pow(math.pow(ds.rg0(),-1/r)-1,-1/k)+l
    @staticmethod
    def pdf(r,k,l,s,x):
        return 1/s*r*k*math.pow((x-l)/s,r*k-1)/math.pow(1+math.pow(x,k),r+1)
    @staticmethod
    def cdf(r,k,l,s,x):
        return math.pow(1+math.pow((x-l)/s,-k),-r)
    @staticmethod
    def median(r,k,l,s):
        return math.pow(math.pow(1/2,-1/r)-1,-1/k)+l
    @staticmethod
    def ppf(r,k,l,s,q):
        return s*math.pow(math.pow(q,-1/r)-1,-1/k)+l
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=burr3.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(1,10),(1,10),(1,10),(1,10)]).x.tolist()
        return {'r':ret[0],'k':ret[1],'l':ret[2],'s':ret[3]}