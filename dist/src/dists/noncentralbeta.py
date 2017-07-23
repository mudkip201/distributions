'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.noncentralchi2.noncentralchi2 as noncentralchi2
import dists.chi2.chi2 as chi2

class noncentralbeta(Distribution):
    @staticmethod
    def random(m,n,mu,sigma):
        return(noncentralchi2.random(mu,sigma,m))/(noncentralchi2.random(mu,sigma,m)+chi2.random(n))