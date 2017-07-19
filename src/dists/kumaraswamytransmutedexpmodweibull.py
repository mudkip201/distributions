'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import numpy as np
from numpy import random as r
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as op

class kumaraswamytransmutedexpmodweibull(Distribution):
    @staticmethod
    def pdf(a,b,g,l,t,aa,bb,x): #note - aa and bb are the Kumaraswamy a and b, and a and b are the alpha and beta
        g=aa*bb*a*(t+g*b*math.pow(x,b-1))*math.exp(-(t*x+g*math.pow(x,b)))*(1+l-2*l*math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a))
        g*=math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a*aa-1)*math.pow(1+l-l*math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a),aa-1)
        g*=math.pow(1-math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a*aa)*math.pow(1+l-l*math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a),aa),bb-1)
        pass
    @staticmethod
    def cdf(a,b,g,l,t,aa,bb,x):
        return 1-math.pow(1-math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a*aa)*math.pow(1+l-l*math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a),aa),bb)
        pass
    @staticmethod
    def random(a,b,g,l,t,aa,bb,epsilon=0.0000001):
        u=r.random()
        diff=1
        x0=1
        while(diff>epsilon):
            xstar=x0-(kumaraswamytransmutedexpmodweibull.cdf(a,b,g,l,t,aa,bb,x0)-u)/kumaraswamytransmutedexpmodweibull.pdf(a,b,g,l,t,aa,bb,x0)
            diff=abs(x0-xstar)
            x0=xstar
        return x0
    @staticmethod
    def median(a,b,g,l,t,aa,bb,epsilon=0.0000001):
        return kumaraswamytransmutedexpmodweibull.ppf(a,b,g,l,t,aa,bb,0.5,epsilon)
    @staticmethod
    def ppf(a,b,g,l,t,aa,bb,q,epsilon=0.0000001):
        diff=1
        x0=1
        while(diff>epsilon):
            xstar=x0-(kumaraswamytransmutedexpmodweibull.cdf(a,b,g,l,t,aa,bb,x0)-q)/kumaraswamytransmutedexpmodweibull.pdf(a,b,g,l,t,aa,bb,x0)
            diff=abs(x0-xstar)
            x0=xstar
        return x0
