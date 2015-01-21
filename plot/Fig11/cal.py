#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc

rc('text', usetex=True)

# read data
#data = np.loadtxt('mass04.txt')
data = np.loadtxt('mass24.txt')

(n,m) = data.shape

i = 0
#for i in np.arange(0,n):
while (i < n):
	ave = 0.0
	h = 0
	for j in np.arange(0,n):
		if (data[i,0] == data[j,0]):
			ave = ave+data[j,1]
			h = h+1
	ave = ave/h
	print data[i,0], ave
	i = i+h

