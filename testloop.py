# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 20:09:04 2016

@author: Ben
"""

from numpy import *


m = meshgrid(range(10), range(10))
xs = m[0]
ys = m[1] 

for ik1,ik2 in zip(xs,ys):
    print ik1,ik2