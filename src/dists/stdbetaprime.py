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
import dists.gamma.gamma as gamma

class stdbetaprime(Distribution): #standard beta prime
    @staticmethod
    def random(a,b):
        c=gamma.random(b,1)
        while(c==0):
            c=gamma.random(b,1)
        return gamma.random(a,1)/c
