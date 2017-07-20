'''
Created on Jul 17, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.kumaraswamytransmutedexpmodweibull.kumaraswamytransmutedexpmodweibull as kumaraswamytransmutedexpmodweibull

class kumaraswamymodweibull(Distribution):
    @staticmethod
    def pdf(a,b,t,g,bb,x):
        return kumaraswamytransmutedexpmodweibull.pdf(1, bb, g, 0, t, a, b, x)
    @staticmethod
    def cdf(a,b,t,g,bb,x):
        return kumaraswamytransmutedexpmodweibull.cdf(1, bb, g, 0, t, a, b, x)
    @staticmethod
    def random(a,b,t,g,bb):
        return kumaraswamytransmutedexpmodweibull.random(1, bb, g, 0, t, a, b)
    @staticmethod
    def median(a,b,t,g,bb):
        return kumaraswamytransmutedexpmodweibull.median(1, bb, g, 0, t, a, b)
    @staticmethod
    def ppf(a,b,t,g,bb):
        return kumaraswamytransmutedexpmodweibull.ppf(1, bb, g, 0, t, a, b)