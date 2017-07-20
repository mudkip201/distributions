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

class moyal(Distribution):
    @staticmethod
    def random():
        x1=ds.rg0()
        x2=ds.rg0()
        y=math.pi*x1-math.pi/2
        h=x2*0.912
        z=math.tan(y)
        hy=1/math.sqrt(2*math.pi)*1/(math.cos(y)^2)*math.exp(-(math.tan(y)+math.exp(-math.tan(y)))/2)
        while(h>hy):
            x1=ds.rg0()
            x2=ds.rg0()
            y=math.pi*x1-math.pi/2
            h=x2*0.912
            z=math.tan(y)
            hy=1/math.sqrt(2*math.pi)*1/(math.cos(y)^2)*math.exp(-(math.tan(y)+math.exp(-math.tan(y)))/2)
        return z
