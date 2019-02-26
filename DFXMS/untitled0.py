# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 11:15:28 2018

@author: Yong
"""

import os

f = open("pFind.spectra", 'r").read
print(f[1].strip().split("\t"))
# for line in f[1:]:
    