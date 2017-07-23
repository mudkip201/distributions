'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class transmutedexpmodweibull(Distribution):
    @staticmethod
    def pdf(a,b,g,l,t,x):
        g=a*(t+g*b*math.pow(x,b-1))*math.exp(-(t*x+g*math.pow(x,b)))
        g*=math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a-1)
        g*=(1+l-2*l*math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a))
        return g
    @staticmethod
    def cdf(a,b,g,l,t,x):
        G=math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a)
        G*=(1+l-l*math.pow(1-math.exp(-(t*x+g*math.pow(x,b))),a))
        return G