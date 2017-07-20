'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''
import dists.Distribution.Distribution as Distribution
import math
import numpy as np
import scipy.optimize as op
import dists.gamma.gamma as gamma

class amoroso(Distribution):
    @staticmethod
    def random(a,theta,aa,bb):
        return a+theta*math.pow(gamma.random(aa,1),1/bb)
    @staticmethod
    def pdf(a,theta,aa,bb,x):
        #print(a,theta,aa,bb,x)
        if(theta>0 and x<a):
            raise ValueError("x must be greater than a")
        if(theta<0 and x>a):
            raise ValueError("x must be less than a")
        aaa=1/math.gamma(aa)
        bbb=abs(bb/theta)
        ccc=math.pow((x-a)/theta,aa*bb-1)
        ddd=math.exp(-math.pow((x-a)/theta,bb))
        return aaa*bbb*ccc*ddd
    @staticmethod
    def mean(a,theta,aa,bb):
        if(aa+1/bb>=0):
            return a+theta*math.gamma(aa+1/bb)/math.gamma(aa)
        return None
    @staticmethod
    def mode(a,theta,aa,bb):
        if(aa*bb>=1):
            return a+theta*math.pow(aa-1/bb,1/bb)
        if(aa*bb<=1):
            return a
    @staticmethod
    def variance(a,theta,aa,bb):
        if(aa+2/bb>=0):
            return theta**2*(math.gamma(aa+2/bb)/math.gamma(aa)-(math.gamma(aa+1/bb)/math.gamma(aa))**2)
        return None
    @staticmethod
    def stddev(a,theta,aa,bb):
        if(aa+2/bb>=0):
            return theta*math.sqrt(math.gamma(aa+2/bb)/math.gamma(aa)-(math.gamma(aa+1/bb)/math.gamma(aa))**2)
        return None
    @staticmethod
    def mle(x): #not working
        args0=[1,0.3,1.05,0.6]
        cons=[]
        for i in x:
            cons.append(lambda y: np.sign(y[1])*i-y[0])
        cons.append(lambda y: min(x)-y[0])
        cons.append(lambda y: y[2])
        def mlefunc(args_):
            tomin=1
            print(args_)
            for i in x:
                tomin+=math.log(amoroso.pdf(args_[0],args_[1],args_[2],args_[3],i))
            return -tomin

        ret=op.slsqp.fmin_slsqp(mlefunc,args0,ieqcons=cons,iprint=2)
        return {'a':ret[0],'theta':ret[1],'aa':ret[2],'bb':ret[3]}