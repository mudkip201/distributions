'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.halfgenpearson7.halfgenpearson7 as halfgenpearson7

class halflaha(Distribution):
    @staticmethod
    def random(a,s):
        return halfgenpearson7.random(a,s,1,4)