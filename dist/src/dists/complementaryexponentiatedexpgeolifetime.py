'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class complementaryexponentiatedexpgeolifetime(Distribution):
    @staticmethod
    def random(a,l,t):
        u=ds.rg0()
        return -(math.log(1-math.pow(u/(t*(1-u)+u),1/a)))/l
    @staticmethod
    def pdf(a,l,t,x):
        return a*l*t*math.exp(-l*x)*math.pow(1-math.exp(-l*x),a-1)/math.pow(1-(1-t)*math.pow(1-math.exp(-l*x),a),2)
    @staticmethod
    def cdf(a,l,t,x):
        z=1-math.exp(-l*x)
        return 1-(1-math.pow(z,a))/(1-(1-t)*math.pow(z,a))
    @staticmethod
    def median(a,l,t):
        return -(math.log(1-math.pow(1/(t+1),1/a)))/l
    @staticmethod
    def ppf(a,l,t,q):
        return -(math.log(1-math.pow(q/(t*(1-q)+q),1/a)))/l
    @staticmethod
    def mle(x):#not in eclipse
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=complementaryexponentiatedexpgeolifetime.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        '''
        ret=op.differential_evolution(mlefunc,[(0.01,10),(0.01,10),(0.01,10)]).x.tolist()
        return {'a':ret[0],'l':ret[1],'t':ret[2]}
'''