'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class transmutedweibulllomax(Distribution):
    @staticmethod
    def random(a,b,aa,bb,l):
        u=ds.rg0()
        D=0
        if(l==0):
            D=u
        else:
            D=((1+l)-math.sqrt((1+l)**2-4*l*u))/(2*l)
        return bb*(math.pow(math.pow(math.log(math.pow(1-D,-1/a)),1/b)+1,1/aa)-1)
    @staticmethod
    def pdf(a,b,aa,bb,l,x):
        return a*b*aa/bb*math.pow(1+x/bb,b*a-1)*math.pow(1-math.pow(1+x/bb,-aa),b-1)*math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))*((1+l)-2*l*(1-math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))))
    @staticmethod
    def cdf(a,b,aa,bb,l,x):
        return (1-math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b)))*((1+l)-l*(1-math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))))
    @staticmethod
    def median(a,b,aa,bb,l):
        D=0
        if(l==0):
            D=1/2
        else:
            D=((1+l)-math.sqrt(1+l**2))/(2*l)
        return bb*(math.pow(math.pow(math.log(math.pow(1-D,-1/a)),1/b)+1,1/aa)-1)
    @staticmethod
    def ppf(a,b,aa,bb,l,q):
        D=0
        if(l==0):
            D=q
        else:
            D=((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return bb*(math.pow(math.pow(math.log(math.pow(1-D,-1/a)),1/b)+1,1/aa)-1)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=transmutedweibulllomax.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'aa':ret[2],'bb':ret[3]}