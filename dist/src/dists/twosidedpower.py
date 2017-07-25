'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class twosidedpower(Distribution):
    @staticmethod
    def pdf(a,b,m,n,x):
        if(a<x and x<=m):
            return n/(b-a)*math.pow((x-a)/(m-a),n-1)
        if(m<=x and x<b):
            return n/(b-a)*math.pow((b-x)/(b-m),n-1)
        return 0
    @staticmethod
    def cdf(a,b,m,n,x):
        if(x<a):
            return 0
        if(a<=x and x<=m):
            return (m-a)/(b-a)*math.pow((x-a)/(m-a),n)
        if(m<=x and x<=b):
            return 1-(b-m)/(b-a)*math.pow((b-x)/(b-m),n)
        if(b<x):
            return 1
    @staticmethod
    def mean(a,b,m,n):
        return (a+(n-1)*m+b)/(n+1)
    @staticmethod
    def mode(a,b,m,n):
        if(n>1):
            return n/(b-a)
        if(0<=n and n<1 and a<m and m<b):
            return (a,b)
    @staticmethod
    def variance(a,b,m,n):
        return (b-a)**2*(n-2*(n-1)*(m-a)/(b-a)*(b-m)/(b-a))/((n+2)*(n+1)**2)
    @staticmethod
    def stddev(a,b,m,n):
        return (b-a)*math.sqrt((n-2*(n-1)*(m-a)/(b-a)*(b-m)/(b-a))/((n+2)))/(n+1)