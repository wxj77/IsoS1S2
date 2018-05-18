# -*- coding: utf-8 -*-
"""
Created on Mon May 14 19:06:24 2018

@author: wei
"""

#readcsv from dev

import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
workdir="/home/wei/"

S1area_step=10 #S1:0.1 S2:10
S1rebin=10;#S1:10 S2:1
S1S2=1 #S1:1 S2:2

if (S1S2==1):
    S1area_step=0.1 #S1:0.1 S2:10
    S1rebin=10;#S1:10 S2:0.1
elif (S1S2==2):
    S1area_step=10 #S1:0.1 S2:10
    S1rebin=1;#S1:10 S2:0.1

livetime = np.zeros(1); 
deadtime = np.zeros(1)
numOfS1s =1000
S1 = np.zeros([1,numOfS1s ])
S1_rb = np.zeros([1,numOfS1s /S1rebin])

livetimes= [] # in second?
deadtimes=[] # in second?
S1s = np.zeros([0, numOfS1s ]);
S1s_rb = np.zeros([0, numOfS1s/S1rebin ]);
S1srate = np.zeros([0, numOfS1s ]);
S1srate_rb = np.zeros([0, numOfS1s/S1rebin ]);
S1sarea = np.linspace(1,numOfS1s*S1area_step+1, numOfS1s,endpoint=False); # S1 start from 1, step=0.01: S2 start from 1 step 10
S1sarea_rb = np.linspace(1,numOfS1s*S1area_step+1, numOfS1s/S1rebin,endpoint=False); # S1 start from 1, step=0.01: S2 start from 1 step 10

with open('/home/wei/Downloads/SlackDownloads/IsoS1S2/total_iso_s'+str(S1S2)+'_spect_real1.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in spamreader:
        lt=float(row[2])*60;
        dt=float(row[3])*60;
        livetimes.append(lt);
        deadtimes.append(lt);
        for jj in range(numOfS1s ):
            S1[0][jj] = float(row[jj+4])/S1area_step
        S1s = np.concatenate((S1s, S1), axis=0);      
        S1srate = np.concatenate((S1srate, S1/lt), axis=0);
        for kk in range(numOfS1s/S1rebin):
            if (kk==0):
                S1_rb[0,kk]=np.mean(S1[0,0:kk*S1rebin+S1rebin]/lt)
            else:                
                S1_rb[0,kk]=np.mean(S1[0,kk*S1rebin:kk*S1rebin+S1rebin]/lt)
        S1srate_rb  = np.concatenate((S1srate_rb, S1_rb), axis=0);

livetimes = np.array(livetimes)

deadtimes = np.array(deadtimes)

S1srate_mean1_temp = np.percentile(S1srate[:,1],90)
S1srate_std1_temp = np.mean(S1srate[:,1])
S1srate_90= np.percentile(S1srate[:,1],90)
S1srate_10= np.percentile(S1srate[:,1],10)
S1srate_95= np.percentile(S1srate[:,1],95)
S1srate_05= np.percentile(S1srate[:,1],5)

S1srate_tr = np.copy(S1srate)
S1srate_rb_tr = np.copy(S1srate_rb)

kkk=0
for kk in range(S1srate.shape[0]):
    if ((S1srate[kk,1]<S1srate_05) | (S1srate[kk,1]>S1srate_95)):
 #       print S1srate[kk,1],kk, kkk
        S1srate_tr = np.delete(S1srate_tr, kk-kkk, 0)
        S1srate_rb_tr = np.delete(S1srate_rb_tr, kk-kkk, 0)
        kkk+=1

S1srate_mean = np.zeros(S1srate.shape[1])
S1srate_median = np.zeros(S1srate.shape[1])



for jj in range(S1srate_mean.shape[0]):
    S1srate_mean[jj] = np.mean(S1srate[:,jj])
    S1srate_median[jj] = np.median(S1srate[:,jj])



S1srate_tr_mean = np.zeros(S1srate_tr.shape[1])
S1srate_tr_median = np.zeros(S1srate_tr.shape[1])
S1srate_rb_tr_mean = np.zeros(S1srate_rb_tr.shape[1])
S1srate_rb_tr_median = np.zeros(S1srate_rb_tr.shape[1])


for jj in range(S1srate_tr_mean.shape[0]):
    S1srate_tr_mean[jj] = np.mean(S1srate_tr[:,jj])
    S1srate_tr_median[jj] = np.median(S1srate_tr[:,jj])

for jj in range(S1srate_rb_tr_mean.shape[0]):
    S1srate_rb_tr_mean[jj] = np.mean(S1srate_rb_tr[:,jj])
    S1srate_rb_tr_median[jj] = np.median(S1srate_rb_tr[:,jj])

def loglinear(x,p0,p1):
    return 10**(p0+p1*np.log10(x))

xx=S1sarea_rb[10:]
yy=S1srate_rb_tr_mean[10:]
best_vals, covar = curve_fit(loglinear, xx, yy, p0=[1,1])
xxf=S1sarea_rb[0:]
yyf = loglinear(xxf,best_vals[0], best_vals[1])

fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
ax.plot(S1sarea, S1srate_mean, label ="mean: no data selection")
ax.plot(S1sarea, S1srate_tr_mean, label ="mean: with data selection")
ax.plot(S1sarea_rb+0.5, S1srate_rb_tr_mean, label ="mean: with data selection(average)")
if (S1S2==1):
    ax.plot(xxf, yyf, label =r"$log_{10}y=(%.3f)+(%.3f)*log_{10}x$"%(best_vals[0], best_vals[1]))
ax.set_xlabel("S"+str(S1S2)+" [phd]")
ax.set_ylabel("isolated S"+str(S1S2)+" rate [Hz/phe]")
ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_ylim(0,0.002)
#ax.text(0.5,0.8, r"$log_{10}y=(%.3f)+(%.3f)*x$"%(best_vals[0], best_vals[1]), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
ax.legend()

plt.savefig(workdir+"isoS"+str(S1S2)+"rate.png")


for jj in range(8):
#    print jj+1.5, S1srate_rb_tr_mean[jj]
    print (S1srate_rb_tr_mean[jj],",")










