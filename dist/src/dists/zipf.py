'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import numpy as np
import scipy.stats as st
import scipy.special as sp

class zipf(Distribution):
    @staticmethod
    def random(rho):
        return st.zipf.rvs(rho)
    @staticmethod
    def pdf(rho,x):
        return math.pow(x,-rho-1)/sp.zeta(rho+1)
    @staticmethod
    def mean(rho):
        if rho>1:
            return sp.zeta(rho)/sp.zeta(rho+1)
        return np.Infinity
    @staticmethod
    def variance(rho):
        if rho>2:
            return sp.zeta(rho-1)/sp.zeta(rho+1)-sp.zeta(rho)**2/sp.zeta(rho+1)**2
        return np.Infinity
    @staticmethod
    def stddev(rho):
        if rho>2:
            return math.sqrt(zipf.variance(rho))
        return np.Infinity
    @staticmethod
    def skewness(rho):
        z0=sp.zeta(rho)
        z1=sp.zeta(rho+1)
        zn1=sp.zeta(rho-1)
        zn2=sp.zeta(rho-2)
        if(rho>3):
            return (2*z0**3/z1**3-3*zn1*z0/z1**2+zn2/z1)/(zn1/z1-z0**2/z1**2)
        return np.Infinity
    @staticmethod
    def kurtosis(rho):
        z0=sp.zeta(rho)
        z1=sp.zeta(rho+1)
        zn1=sp.zeta(rho-1)
        zn2=sp.zeta(rho-2)
        zn3=sp.zeta(rho-3)
        if(rho>4):
            return (-3*z0**4+6*zn1*z1*z0**2-4*zn2*z1**2*z0+zn3*z1**3)/(z0**2-zn1*z1)**2
        return np.Infinity