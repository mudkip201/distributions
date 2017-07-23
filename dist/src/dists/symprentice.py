'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.prentice.prentice as prentice

class symprentice(Distribution): #symmetric prentice
    @staticmethod
    def random(l,a):
        return prentice.random(0,l,a,a)