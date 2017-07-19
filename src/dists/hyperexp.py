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
import dists.exponential.exponential as exp

class hyperexp(Distribution): #hyperexponential
    @staticmethod
    def random(lmbdas,probs):
        a=r.random()
        probs2=probs/sum(probs)
        counter=0
        while(a>sum(probs2[:counter+1])):
            counter+=1
        return exp.random(lmbdas[counter])