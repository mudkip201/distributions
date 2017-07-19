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
import dists.genbetaprime.genbetaprime as genbetaprime


class halfgenpearson7(Distribution): #half generalized pearson VII
    @staticmethod
    def random(a,s,bb,m):
        return genbetaprime.random(a,s,1/bb,m-1/bb,bb)
