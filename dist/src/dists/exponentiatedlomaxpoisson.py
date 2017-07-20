'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class exponentiatedlomaxpoisson(Distribution):
    @staticmethod
    def random(a,b,g,l):#bad
        u=ds.rg0()
        return (math.pow(-math.pow(-math.log(-u*(math.exp(l)-1)+math.exp(l))/l,1/a)+1,-1/g)-1)/b
    @staticmethod
    def pdf(a,b,g,l,x): #bad?
        return l*a*g*b*math.exp(l)*math.pow(1-math.pow(1+b*x,-g),a-1)*math.exp(-l*math.pow(1-math.pow(1+b*x,-g),a))/((math.exp(l)-1)*math.pow(1+b*x,g+1))
    @staticmethod
    def cdf(a,b,g,l,x): #bad?
        return math.exp(l)*(1-math.exp(-l*math.pow(1-math.pow(1+b*x,-g),a)))/(math.exp(l)-1)
    @staticmethod
    def median(a,b,g,l):#bad?
        return (math.pow(-math.pow(-math.log((math.exp(l)/2-1)+math.exp(l))/l,1/a)+1,-1/g)-1)/b
    @staticmethod
    def ppf(a,b,g,l,q):#bad?
        return (math.pow(-math.pow(-math.log(-q*(math.exp(l)-1)+math.exp(l))/l,1/a)+1,-1/g)-1)/b
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=exponentiatedlomaxpoisson.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'g':ret[2],'l':ret[3]}
