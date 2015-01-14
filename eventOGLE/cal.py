#!/usr/bin/env python

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc

ogle0 = np.loadtxt('ogleUse.txt')
rate = np.loadtxt('map_event_beta-1.txt')

n1,m1 = ogle0.shape
n2,m2 = rate.shape
#print n1,m1
#print n2,m2

ogle = (np.pi/180.0)*ogle0
#print 180.0*ogle[0,1]/ogle0[0,1]

A = (np.pi/180.0)*0.1
eventRate = 0.0
n = 0
for i in np.arange(0,n1):
#	print ogle[i,0]
		for j in np.arange(0,n2):
				if (np.fabs(rate[j,0]-ogle[i,0]) <= A and np.fabs(rate[j,1]-ogle[i,1]) <= A):
					eventRate = eventRate+rate[j,2]
					n = n+1

print n, eventRate, eventRate/n 

n_star = 340e6

print n_star*eventRate/n
p_psr = (eventRate/n)*120e3/1e9

nu = p_psr*n_star
nu = nu*9

p = 0
n_tot = 0

for n in np.arange(1,20):
	pn = (nu**n)*math.exp(-nu)/math.factorial(n)
	n_tot = n_tot+n*pn
	p = p+pn

print p*100, n_tot
