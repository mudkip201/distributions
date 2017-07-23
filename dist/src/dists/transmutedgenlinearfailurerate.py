'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class transmutedgenlinearfailurerate(Distribution): #transmuted generalized linear failure rate
    @staticmethod
    def random(a,g,l,t):
        i=((1+l)-math.sqrt((1+l)**2-4*l*ds.rg0()))/(2*l)
        return (-t+math.sqrt(t**2-4*g*math.log(1-math.pow(i,1/a))))/(2*g)
    @staticmethod
    def median(a,g,l,t):
        i=((1+l)-math.sqrt((1+l)**2-2*l))/(2*l)
        return (-t+math.sqrt(t**2-4*g*math.log(1-math.pow(i,1/a))))/(2*g)
    @staticmethod
    def ppf(a,g,l,t,q):
        i=((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return (-t+math.sqrt(t**2-4*g*math.log(1-math.pow(i,1/a))))/(2*g)