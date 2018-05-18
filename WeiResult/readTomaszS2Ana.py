# -*- coding: utf-8 -*-
"""
Created on Tue May 15 13:27:08 2018

@author: wei
"""

# this code is to generate the s2 spectrum from Tomasz's analysis

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
workdir="/home/wei/"

S1S2=2
a = np.loadtxt("isolatedS2s_livetime_areas_allS2QualityCuts-fidCuts_180513.txt")
lt = a[0]
S2s =a[1:]

bins=10**(np.linspace(2,4.6, 1000))
value,bins = np.histogram(S2s, bins=bins)
value=value/lt
binsx= 0.5*(bins[1:]+bins[:-1])

def loglinear(x,p0,p1):
    return 10**(p0+p1*np.log10(x))

xx=binsx[100:900]
yy=(value/(bins[1:]-bins[:-1]))[100:900]
best_vals, covar = curve_fit(loglinear, xx, yy, p0=[-7,-1])
xxf=binsx[0:]
yyf = loglinear(xxf,best_vals[0], best_vals[1])

fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)

ax.step(bins[1:], value/(bins[1:]-bins[:-1]), where='pre', color='b', label="isolated S2")
ax.step(bins[:-1], value/(bins[1:]-bins[:-1]), where='post', color='b')
if (S1S2==2):
    ax.plot(xxf, yyf, label =r"$log_{10}y=(%.3f)+(%.3f)*log_{10}x$"%(best_vals[0], best_vals[1]))
ax.set_xlabel("S"+str(S1S2)+" [phd]")
ax.set_ylabel("isolated S"+str(S1S2)+" rate [Hz/phe]")
ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_ylim(0,0.002)
#ax.text(0.5,0.8, r"$log_{10}y=(%.3f)+(%.3f)*x$"%(best_vals[0], best_vals[1]), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
ax.legend()

plt.savefig(workdir+"isoS"+str(S1S2)+"rate.png")


