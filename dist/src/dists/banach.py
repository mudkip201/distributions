'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Dist
import dists.Distribution as ds
import math


class banach(Dist):
    @staticmethod
    def random(x,n):
        p=ds.rg0()
        sum_=0
        counter=0
        while(sum_<n):
            sum_+=math.pow(2,x-2*p)*(math.gamma(2*p-x+1)/(math.gamma(p+1)*math.gamma(p-x+1)))
            counter+=1
        return counter