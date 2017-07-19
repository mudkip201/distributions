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
import dists.halfgenpearson7.halfgenpearson7 as halfgenpearson7

class genpearson7(Distribution): #generalized pearson VII
    @staticmethod
    def random(a,s,m,b):
        n=r.random()
        if(n<0.5):
            return -halfgenpearson7.random(a,s,m,b)
        return halfgenpearson7.random(a,s,m,b)
