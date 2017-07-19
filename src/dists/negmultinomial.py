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

class negmultinomial(Distribution): #negative multinomial
    @staticmethod
    def random(probs, n):
        amt=[]
        for _ in range(len(probs)):
            amt.append(0)
        probs2=probs/sum(probs)
        while(amt[0]<n):
            a=r.random()
            i=0
            while(a>sum(probs2[:i+1])):
                i+=1
            amt[i]+=1
        return amt
