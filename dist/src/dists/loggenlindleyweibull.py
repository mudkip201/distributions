'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op
import dists.gamma.gamma as gamma

class loggenlindleyweibull(Distribution): #log generalized lindley-weibull
    @staticmethod
    def random(a,b,t,g,c):
        u=ds.rg0()
        v1=gamma.random(a,t)
        v2=gamma.random(a+1,t)
        if(u<=t/(b+t)):
            return g*math.pow(v1,1/c)
        return g*math.pow(v2,1/c)
    @staticmethod
    def pdf(a,b,t,g,c,x):
        return c*math.pow(t,a+1)/(g*(b+t)*math.gamma(a+1))*math.pow(x/g,c*a-1)*(a+b*math.pow(x/g,c))*math.exp(-t*math.pow(x/g,c))
    @staticmethod
    def mean(a,b,t,g,c):
        return g*math.pow(t,1-1/c)/((b+t)*math.gamma(a+1))*math.gamma(a+1/c)*(a+b/t*(a+1/c))
    @staticmethod
    def mode(a,b,t,g,c):
        return g*math.pow((c*(b-t)*a+b*(c-1))/(2*c*t*b)+math.sqrt(math.pow(-c*(b-t)*a-b*(c-1),2)-4*b*c*t*(-c*a*a+a))/(2*c*t*b),1/c)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=loggenlindleyweibull.pdf(args_[0],args_[1],args_[2],args_[3],args_[4],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'t':ret[2],'g':ret[3],'c':ret[4]}