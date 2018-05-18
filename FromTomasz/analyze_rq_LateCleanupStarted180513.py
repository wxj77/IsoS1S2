
#import read_root as rr
import numpy as n
import os, sys, glob
import re
#import fitgausshist as fg
import cPickle as cp
import time
#import get_html_data as gethtml
import gc

import hist2d

from scipy.interpolate import interp1d
import extrap1d

#from scipy.optimize import leastsq
#from scipy import stats

import matplotlib.pyplot as pyp
import matplotlib as mlib		# for color scale
from matplotlib import mpl		# for colorbar
from matplotlib.ticker import Formatter	# to convert xaxis label to dates

import fitgausshist as fg

#from datetime import datetime

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Function definitions
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

def get_isolated(which, rq, inds, s2cutoff=55, dos2100=False):
	"""
	Find isolated S1s or S2s. 
	"""
	rqnum = len(rq)
	# figure out which pulse class we're doing
	if which == 's1':
		class1 = 1
		class2 = 2
	elif which == 's2':
		class1 = 2
		class2 = 1
	else:
		print "Only 's1' and 's2' are allowed keywords for 'which'"
		raise(ValueError)
	#
	inds_iso = []				# without S1s or S2s
	inds_iso_unknown = []		# without S1s or S2s or unknowns
	if dos2100: inds_iso_100phe = []		# without S1s or S2s. Use 100 phe as a cutoff
	for i in range(rqnum):
		tmprow, tmpcol = [], []
		tmprowunk, tmpcolunk = [], []
		if dos2100: tmprow100, tmpcol100 = [], []
		for j in range(inds[i][0].size):
			temps2flags = rq[i]['pulse_classification'][inds[i][0][j]] == class1
			temps1flags = rq[i]['pulse_classification'][inds[i][0][j]] == class2
			if which == 's1':
				temps1flags = temps1flags & (rq[i]['pulse_area_phe'][inds[i][0][j]] > s2cutoff)
			elif which == 's2':
				temps2flags = temps2flags & (rq[i]['pulse_area_phe'][inds[i][0][j]] > s2cutoff)
			if temps2flags.sum() == 1 and temps1flags.sum() == 0:
				# found isolated S2 (not accounting for unkowns)
				tmprow.append(inds[i][0][j])
				tmpcol.append(inds[i][1][j])
				if (rq[i]['pulse_classification'][inds[i][0][j]] == 5).sum() == 0:
					# also isolated from unkowns
					tmprowunk.append(inds[i][0][j])
					tmpcolunk.append(inds[i][1][j])
			if dos2100: 
				# use a 100 phe cutoff for S2 definition to determine the isolated S2 rate
				temp100flags = rq[i]['pulse_area_phe'][inds[i][0][j]] > 100
				if (temps2flags & temp100flags).sum() == 1 and \
					temps1flags.sum() == 0:# and \
					#(rq[i]['pulse_classification'][inds[i][0][j]] == 5).sum() == 0:
					tmprow100.append(inds[i][0][j])
					tmpcol100.append(inds[i][1][j])		
		inds_iso.append(( n.array(tmprow), n.array(tmpcol) ))
		inds_iso_unknown.append(( n.array(tmprowunk), n.array(tmpcolunk) ))
		if dos2100: inds_iso_100phe.append(( n.array(tmprow100), n.array(tmpcol100) ))
	# return indecies
	if dos2100:
		return inds_iso, inds_iso_unknown, inds_iso_100phe
	else:
		return inds_iso, inds_iso_unknown


def get_golden_inds(rq, rqnum, s2onlyflags=None, secut=None, second_smaller_s2cut=None, selected_s1_s2_key = 'selected_s1_s2', allows1num = 1):
	"""
	Return indecies of golden events
	
	S2s are defined only if larger than secut. 
	Second S2 is allowed as long as it's smaller than the first and smaller than
	second_smaller_s2cut.
	
	"""
	fidgoldflags_s1 = []
	fidgoldflags_s2 = []
	fidgoldinds_s1 = []
	fidgoldinds_s2 = []
	for i in range(rqnum):
		if s2onlyflags == None:
			ts2onlyflags = rq[i]['pulse_classification'] == 2
		else:
			ts2onlyflags = s2onlyflags[i]
		# Check if S2s are required to be larger than secut
		if secut == None:
			# Use S2s as they appear. 
			senonflag = n.zeros(rq[i]['pulse_classification'].shape, dtype=n.bool)+True
		else:
			# flag pulses larger than secut
			senonflag = rq[i]['pulse_area_phe'] > secut
		sel12flag = rq[i][selected_s1_s2_key] == 1
		# mark events with paired S1s and S2s where the S2s are larger secut
		temprowstemp, tempcols = n.where(ts2onlyflags & sel12flag & senonflag)
		temprowstemp = n.unique(temprowstemp)
		# for some reason, the golden definition for the second Kr set is a complete
		# failure. There are a buch of multiple scatters (multiple S2s) in hear. Deal
		# with it.
		temps1col, temps2col, temprows = [], [], []
		for j in range(temprowstemp.size): 
			# check the number of s2s
			s2flagsperevt = (rq[i]['pulse_classification'][temprowstemp[j]] == 2) & senonflag[temprowstemp[j]]
			s2flagsperevtsum = s2flagsperevt.sum()
			if second_smaller_s2cut != None and s2flagsperevtsum > 1:
				# More than 1 S2 was found but we are going to keep this event
				# if the additinal S2 are small enough.
				multis2inds = n.where(s2flagsperevt)[0]
				if n.all( rq[i]['pulse_area_phe'][temprowstemp[j]][multis2inds[1:]] < rq[i]['pulse_area_phe'][temprowstemp[j]][multis2inds[0]] ) and  n.all( rq[i]['pulse_area_phe'][temprowstemp[j]][multis2inds[1:]] < second_smaller_s2cut ):
					s2flagsperevtsum = 1
			# check number of S1ss
			s1flagsperevt = rq[i]['pulse_classification'][temprowstemp[j]] == 1
			nums1flgasperevt = s1flagsperevt.sum()
			if nums1flgasperevt > 0 and nums1flgasperevt <= allows1num and s2flagsperevtsum == 1: 
				# If we got here, we have golden events.
				temprows.append(temprowstemp[j])
				temps1col.append(n.where(s1flagsperevt)[0][0])
				temps2col.append(n.where(s2flagsperevt)[0][0])
		temprows = n.array(temprows)
		fidgoldinds_s1.append((temprows, n.array(temps1col)))
		fidgoldinds_s2.append((temprows, n.array(temps2col)))
		#
		fidgoldflags_s1.append(n.zeros(rq[i]['pulse_area_phe'].shape, dtype=n.bool))
		if fidgoldinds_s1[i][0].size > 0: fidgoldflags_s1[i][fidgoldinds_s1[i]] = True
		fidgoldflags_s2.append(n.zeros(rq[i]['pulse_area_phe'].shape, dtype=n.bool))
		if fidgoldinds_s2[i][0].size > 0: fidgoldflags_s2[i][fidgoldinds_s2[i]] = True
	return fidgoldflags_s1, fidgoldflags_s2, fidgoldinds_s1, fidgoldinds_s2


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Load the RQs

#pickledir = '/u/gl/tomaszbi/lux/studies/160623_run4coincidenceBack/myanalysis'
pickledir = '/home/tomaszbi/lux/studies/160623_run4coincidenceBack/myanalysis/processed_rq/large'
dsets = ['lux10_201411', 'lux10_201507', 'lux10_201604']
rqnum = len(dsets)


# smaller sets for prototyping - first 20 files from each data set (40 data sets per collection)
#oldtimestamps = ['20160628T230527', '20160628T230527', '20160628T230527']

# Large set
oldtimestamps = ['20160629T004205', '20160629T004236', '20160629T004308']

rq, evtnum, oevtnum, evtdifftime = [], [], [], []
for i in range(rqnum):
	rq.append({})
	filename = os.path.join(pickledir, '%s_backgroundRQs_%s_*.npy' %(dsets[i], oldtimestamps[i]))
	filelist = glob.glob(filename)
	for sfile in filelist:
		key = os.path.basename(sfile)[43:-4]
		rq[-1][key] = n.load(sfile)
	evtnum.append(rq[-1]['evtnum-oevtnum-evtdifftime'][0])
	oevtnum.append(rq[-1]['evtnum-oevtnum-evtdifftime'][1])
	evtdifftime.append(rq[-1]['evtnum-oevtnum-evtdifftime'][2])

evtnum, oevtnum, evtdifftime = n.array(evtnum), n.array(oevtnum), n.array(evtdifftime)

# convert x and y to r and phi
r_cm_sq = [ rq[i]['x_cm']**2. + rq[i]['y_cm']**2. for i in range(rqnum)]
r_cm = [ n.sqrt(r_cm_sq[i]) for i in range(rqnum)]
phi = [n.where(rq[i]['y_cm'] >=0, n.arccos(rq[i]['x_cm']/r_cm[i]), -n.arccos(rq[i]['x_cm']/r_cm[i])) for i in range(rqnum)]


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Load Dev's S2 spectrum

devbes2left, devbes2right, devbins2rate, devbins2percenterr = n.loadtxt('/home/tomaszbi/lux/studies/160623_run4coincidenceBack/20160623_s2_rate.csv')


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# The analysis begins

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# compute the livetimes
livetimes = []
evtnums = []
for i in range(rqnum):
	livetimes.append((rq[i]['livetime_end_samples'] - rq[i]['livetime_latch_samples']).sum()/1e8)
	evtnums.append(rq[i]['pulse_classification'].shape[0])

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# select s2

s2flags = []
s2inds = []
for i in range(rqnum):
	s2flags.append(rq[i]['pulse_classification'] == 2)
	s2inds.append(n.where(s2flags[i]))

s1flags = []
s1inds = []
for i in range(rqnum):
	s1flags.append(rq[i]['pulse_classification'] == 1)
	s1inds.append(n.where(s1flags[i]))


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# find isolated S2s
s2inds_iso, s2inds_iso_unknown, s2inds_iso_100phe = get_isolated('s2', rq, s2inds, dos2100=True)

# find isolated S1s
s1inds_iso, s1inds_iso_unknown = get_isolated('s1', rq, s1inds)

print s2inds_iso[0][0].size, s2inds_iso[1][0].size, s2inds_iso[2][0].size
print s2inds_iso_unknown[0][0].size, s2inds_iso_unknown[1][0].size
print s2inds_iso_100phe[0][0].size, s2inds_iso_100phe[1][0].size


print s2inds_iso[0][0].size/livetimes[0], s2inds_iso[1][0].size/livetimes[1], s2inds_iso[2][0].size/livetimes[2]

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Apply the bad area cut


#For Good Area <= 253 phd: Bad Area < 80 phd
#For Good Area > 253 phd: Bad Area < 80 * (((10^(-2.4)) * goodarea)**0.4668)


# For S1s it's simple since they're always smaller than 253. Of course, when
# paired with some random S2 the good area may get larger than 253 in which case
# the acceptable bad area gets larger. I will later do a full sim where I input
# an S1 spectrum and generate random pairs of S1s and S2s.

# Do for S1
passed_s1_iso_bad_area_flags = []
for i in range(rqnum):
	good_area_s1_iso = rq[i]['pulse_area_phe'][s1inds_iso[i]]
	bad_area_s1_iso = rq[i]['full_evt_area_phe'][s1inds_iso[i][0]] - good_area_s1_iso
	passed_s1_iso_bad_area_flags.append(n.zeros(s1inds_iso[i][0].size, dtype=n.bool))
	passed_s1_iso_bad_area_flags[i][(good_area_s1_iso <= 253) & (bad_area_s1_iso < 80)] = True
	passed_s1_iso_bad_area_flags[i][(good_area_s1_iso > 253) & (bad_area_s1_iso < 80  * (((10**(-2.4)) * good_area_s1_iso)**0.4668))] = True

# Do for S2s
passed_s2_iso_bad_area_flags = []
for i in range(rqnum):
	good_area_s2_iso = rq[i]['pulse_area_phe'][s2inds_iso[i]]
	bad_area_s2_iso = rq[i]['full_evt_area_phe'][s2inds_iso[i][0]] - good_area_s2_iso
	passed_s2_iso_bad_area_flags.append(n.zeros(s2inds_iso[i][0].size, dtype=n.bool))
	passed_s2_iso_bad_area_flags[i][(good_area_s2_iso <= 253) & (bad_area_s2_iso < 80)] = True
	passed_s2_iso_bad_area_flags[i][(good_area_s2_iso > 253) & (bad_area_s2_iso < 80  * (((10**(-2.4)) * good_area_s2_iso)**0.4668))] = True


# print the effect of the bad area cut
for i in range(rqnum): print float(passed_s1_iso_bad_area_flags[i].sum()) / s1inds_iso[i][0].size


for i in range(rqnum): print float(passed_s2_iso_bad_area_flags[i].sum()) / s2inds_iso[i][0].size


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# S2 histogram after bad area cut
bes = 10**n.linspace(2,4,10)
bes = n.arange(100,1e4,10)
i=0


a=pyp.loglog()
alldata = []
lss = ['-', '--', ':']

for i in range(rqnum):
	labeltxt = 'Isolated S2' if i == 0 else ''
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, 'b', drawstyle='steps-pre', label=labeltxt, ls=lss[i])
	a=pyp.plot(bes[:-1], histdata, 'b', drawstyle='steps-post', ls=lss[i])

for i in range(rqnum):
	labeltxt = 'Isolated S2 + Bad Area Cut' if i == 0 else ''
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, 'g', drawstyle='steps-pre', label=labeltxt, ls=lss[i])
	a=pyp.plot(bes[:-1], histdata, 'g', drawstyle='steps-post', ls=lss[i])

a=pyp.plot(devbes2right, devbins2rate, 'r', drawstyle='steps-pre', label="S2 spectrum From DATs")
a=pyp.plot(devbes2left, devbins2rate, 'r', drawstyle='steps-post')
a=pyp.xlabel('S2 (phd)', fontsize=18)
a=pyp.ylabel('Rate/Bin (Hz/10phd)', fontsize=18)
a=pyp.legend()
pyp.show()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# S1 histogram
bes = n.arange(1,100,1)
i=0

histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s1inds_iso[i]], bins=bes)
histdata = n.float64(histdata) / livetimes[i]


a=pyp.loglog()
a=pyp.plot(bes[1:], histdata, drawstyle='steps-pre')
a=pyp.plot(bes[:-1], histdata, drawstyle='steps-post')
pyp.show()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# do a quick check of spike count vs area

# ----------- Apply a correction to spike counts - 2D
# read in the correction to spike counts
tempmat = n.loadtxt('/home/tomaszbi/lux/studies/140716_spikeCout_vs_pulseArea/model_pileup/spike_count_correction_v10.dat')
corr_spike_count = tempmat[:,0]
corr_p1s = tempmat[:,1:]
corr_dict = {}
corr_min_spike = corr_spike_count.min()
corr_max_spike = corr_spike_count.max()
corr_default_pars = n.array([1.,0,0,0])	# apply no correction
for j in range(corr_spike_count.size):
	corr_dict[int(corr_spike_count[j])] = corr_p1s[j]

def get_corr_params(corr_dict, spike):
	if spike < corr_min_spike:
		return corr_dict[corr_min_spike]	# smallest correction
	if spike > corr_max_spike:
		return corr_dict[corr_max_spike]
	return corr_dict[int(spike)]

def ratiomodelfunc(p,area):
	return p[0] + p[1]*area + p[2]*area**2. + p[3]*area**3.

def apply_corrections(spikes, areas):
	retarr = n.zeros(spikes.size, dtype=n.float32)
	for j in range(spikes.size):
		retarr[j] = spikes[j] * ratiomodelfunc(get_corr_params(corr_dict, spikes[j]), areas[j])
	return retarr

# apply the correction
pulse_cspike_phe_bot = []
pulse_cspike_phe_top = []
pulse_cspike_phe = []
for i in range(rqnum):
	# correct the spike count for overlap bias
	pulse_cspike_phe_bot.append(n.empty(passed_s1_iso_bad_area_flags[i].sum(), dtype=n.float32) + n.nan)
	pulse_cspike_phe_top.append(n.empty(passed_s1_iso_bad_area_flags[i].sum(), dtype=n.float32) + n.nan)
	# select correctable spikes
	useinds = n.where(
		(rq[i]['spike_count_bot'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]] <= corr_max_spike) &
		(rq[i]['spike_count_top'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]] <= corr_max_spike) & (rq[i]['pulse_area_phe'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]] <= 130))
	pulse_cspike_phe_bot[i][useinds] = apply_corrections(rq[i]['spike_count_bot'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]][useinds], rq[i]['pulse_area_phe_bot'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]][useinds])
	pulse_cspike_phe_top[i][useinds] = apply_corrections(rq[i]['spike_count_top'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]][useinds], rq[i]['pulse_area_phe_top'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]][useinds])
	# total
	pulse_cspike_phe.append(pulse_cspike_phe_bot[i] + pulse_cspike_phe_top[i])

# -----------

i = 0
fininds = n.where(n.isfinite(pulse_cspike_phe[i]))
a=pyp.plot(rq[i]['pulse_area_phe'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]][fininds],pulse_cspike_phe[i][fininds], '.')
a=pyp.xlim([-.2,10])
a=pyp.ylim([-.2,10])
pyp.show()


i = 0
fininds = n.where(n.isfinite(pulse_cspike_phe[i]))
a=pyp.plot(pulse_cspike_phe[i][fininds], rq[i]['pulse_area_phe'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]][fininds], '.')
a=pyp.xlim([-.2,10])
a=pyp.ylim([-.2,10])
pyp.show()


i = 0
fininds = n.where(n.isfinite(pulse_cspike_phe[i]))
a=pyp.plot(pulse_cspike_phe[i][fininds], rq[i]['pulse_area_phe'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]][fininds]/pulse_cspike_phe[i][fininds], '.')
a=pyp.xlim([-.2,10])
a=pyp.ylim([-.2,10])
pyp.show()




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Apply Evan's cuts to see how efficient they are


sigThr_flags = []
aftSigRatio_flags = []
chi2cut_flags = []
for i in range(rqnum):
	# Sigma threshold
	sigThr_flags.append(rq[i]['gaus_fit_sigma_samples'][s2inds_iso[i]] > 35)
	# aft-sigma ratio
	aftSigRatio = (rq[i]['aft_t1_samples'][s2inds_iso[i]] - rq[i]['aft_t0_samples'][s2inds_iso[i]]) / rq[i]['gaus_fit_sigma_samples'][s2inds_iso[i]]
	flags_below = \
		( (rq[i]['pulse_area_phe'][s2inds_iso[i]] < 1200) &
			(aftSigRatio < (2.4+rq[i]['pulse_area_phe'][s2inds_iso[i]]*((3.239212144-2.4)/1200.0))) ) | \
		( (rq[i]['pulse_area_phe'][s2inds_iso[i]] >= 1200) &
			(aftSigRatio < (3.5383-rq[i]['pulse_area_phe'][s2inds_iso[i]]*2.69535e-4 + rq[i]['pulse_area_phe'][s2inds_iso[i]]**2.0 * 1.69126e-08) ) )
	flags_above = aftSigRatio > (1.2+(0.3/10000.0)*rq[i]['pulse_area_phe'][s2inds_iso[i]])
	aftSigRatio_flags.append(flags_below & flags_above)
	# position fit
	chi2cut_flags.append(( (rq[i]['pulse_area_phe'][s2inds_iso[i]] <= 4200) & (rq[i]['chi2'][s2inds_iso[i]]<(50.0+((100.0/4200.0)*rq[i]['pulse_area_phe'][s2inds_iso[i]])))) | \
		( (rq[i]['pulse_area_phe'][s2inds_iso[i]] > 4200) & (rq[i]['chi2'][s2inds_iso[i]]<300) ))


total_s2_flags = [passed_s2_iso_bad_area_flags[i] & sigThr_flags[i] & aftSigRatio_flags[i] & chi2cut_flags[i] for i in range(rqnum)]

i=2
# demonstrate the aftSigRatio_flags
a=pyp.plot(rq[i]['pulse_area_phe'][s2inds_iso[i]], aftSigRatio, '.')
a=pyp.plot(rq[i]['pulse_area_phe'][s2inds_iso[i]][aftSigRatio_flags[i]], aftSigRatio[aftSigRatio_flags[i]], 'x')
a=pyp.ylim([0,8])
a=pyp.xlim([0,10000])
pyp.show()

i=0
# demonstrate the chi2cut_flags
a=pyp.plot(rq[i]['pulse_area_phe'][s2inds_iso[i]], rq[i]['chi2'][s2inds_iso[i]], '.')
a=pyp.plot(rq[i]['pulse_area_phe'][s2inds_iso[i]][chi2cut_flags[i]], rq[i]['chi2'][s2inds_iso[i]][chi2cut_flags[i]], 'x')
a=pyp.ylim([0,800])
a=pyp.xlim([0,10000])
pyp.show()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# S2 histogram after bad area cut and evan's cuts
"""
bes = 10**n.linspace(2,4,10)
bes = n.arange(100,1e4,10)
i=2


a=pyp.loglog()
histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]], bins=bes)
histdata = n.float64(histdata) / livetimes[i]
a=pyp.plot(bes[1:], histdata, 'b', lw=2, drawstyle='steps-pre', label='Isolated S2')
a=pyp.plot(bes[:-1], histdata, 'b', lw=2, drawstyle='steps-post')
histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i]], bins=bes)
histdata = n.float64(histdata) / livetimes[i]
a=pyp.plot(bes[1:], histdata, 'g', lw=2, drawstyle='steps-pre', label='Isolated S2 + Bad Area Cut')
a=pyp.plot(bes[:-1], histdata, 'g', lw=2, drawstyle='steps-post')

histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & sigThr_flags[i]], bins=bes)
histdata = n.float64(histdata) / livetimes[i]
a=pyp.plot(bes[1:], histdata, '--', color='c', drawstyle='steps-pre', label='Isolated S2 + Bad Area Cut + Sigma Thrsh.')
a=pyp.plot(bes[:-1], histdata, '--', color='c',  drawstyle='steps-post')

histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & aftSigRatio_flags[i]], bins=bes)
histdata = n.float64(histdata) / livetimes[i]
a=pyp.plot(bes[1:], histdata, '--', color='orange', drawstyle='steps-pre', label='Isolated S2 + Bad Area Cut + AFT-Sigma Ratio')
a=pyp.plot(bes[:-1], histdata, '--', color='orange',  drawstyle='steps-post')

histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & chi2cut_flags[i]], bins=bes)
histdata = n.float64(histdata) / livetimes[i]
a=pyp.plot(bes[1:], histdata, '--', color='m', drawstyle='steps-pre', label=r'Isolated S2 + Bad Area Cut + $\chi^{2}$')
a=pyp.plot(bes[:-1], histdata, '--', color='m',  drawstyle='steps-post')


histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & sigThr_flags[i] & aftSigRatio_flags[i] & chi2cut_flags[i]], bins=bes)
histdata = n.float64(histdata) / livetimes[i]
a=pyp.plot(bes[1:], histdata, '-', color='k', lw=1.5, drawstyle='steps-pre', label=r'Isolated S2 + Bad Area Cut + All S2 Cuts')
a=pyp.plot(bes[:-1], histdata, '-', color='k', lw=1.5,  drawstyle='steps-post')


a=pyp.plot(devbes2right, devbins2rate, 'r:', lw=1.5, drawstyle='steps-pre', label="Dev's S2 spectrum")
a=pyp.plot(devbes2left, devbins2rate, 'r:', lw=1.5, drawstyle='steps-post')
a=pyp.legend()

a=pyp.ylabel('Rate/Bin (Hz/10phd)', fontsize=18)
a=pyp.xlabel('Raw S2 Area (phd)', fontsize=18)
pyp.show()
"""
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# S2 histogram after bad area cut and evan's cuts - Nicer version
bes = 10**n.linspace(2,4,40)

a=pyp.figure(1, figsize=(10,8))
for i in range(rqnum):
	a=pyp.suptitle('Isolated S2 Distributions', fontsize=18)
	a=pyp.subplot(2,2,i+1)
	a=pyp.grid(which='both')
	a=pyp.loglog()
	a=pyp.title(dsets[i])
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, 'b', lw=2, drawstyle='steps-pre', label='All')
	a=pyp.plot(bes[:-1], histdata, 'b', lw=2, drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, 'g', lw=2, drawstyle='steps-pre', label='Bad Area Cut')
	a=pyp.plot(bes[:-1], histdata, 'g', lw=2, drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & sigThr_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, '--', color='c', drawstyle='steps-pre', label='Bad Area Cut + Sigma Thrsh.')
	a=pyp.plot(bes[:-1], histdata, '--', color='c',  drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & aftSigRatio_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, '--', color='orange', drawstyle='steps-pre', label='Bad Area Cut + AFT-Sigma Ratio')
	a=pyp.plot(bes[:-1], histdata, '--', color='orange',  drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & chi2cut_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, '--', color='m', drawstyle='steps-pre', label=r'Bad Area Cut + $\chi^{2}$')
	a=pyp.plot(bes[:-1], histdata, '--', color='m',  drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & sigThr_flags[i] & aftSigRatio_flags[i] & chi2cut_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, '-', color='k', lw=1.5, drawstyle='steps-pre', label=r'Bad Area Cut + All S2 Cuts')
	a=pyp.plot(bes[:-1], histdata, '-', color='k', lw=1.5,  drawstyle='steps-post')
	#
	if i == 2:
		a=pyp.ylabel('Rate/Bin (Hz)', fontsize=18, y=1.1)
		a=pyp.xlabel('Raw S2 Area (phd)', fontsize=18, x=1.1)
	if i == 2:
		a=pyp.legend(loc = 'upper right', bbox_to_anchor = (1.9, 1.1))

pyp.show()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# S2 radial histogram after bad area cut and evan's cuts
areabes = 10**n.linspace(2,4,5)
radsqbes = n.linspace(0,30**2.,10)

a=pyp.figure(1, figsize=(10,8))
for i in range(rqnum):
	a=pyp.suptitle('Isolated S2 Distributions', fontsize=18)
	a=pyp.subplot(2,2,i+1)
	a=pyp.grid(which='both')
	a=pyp.loglog()
	a=pyp.title(dsets[i])
	tempflags = (rq[i]['pulse_area_phe'][s2inds_iso[i]] >= areabes[0]) & (rq[i]['pulse_area_phe'][s2inds_iso[i]] < areabes[1])
	histdata, nothing = n.histogram(r_cm_sq[i][s2inds_iso[i]][tempflags], bins=radsqbes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(radsqbes[1:], histdata, 'b', lw=2, drawstyle='steps-pre', label='%.1e to %.1e' %(areabes[0], areabes[1]))
	a=pyp.plot(radsqbes[:-1], histdata, 'b', lw=2, drawstyle='steps-post')
	
		#
	if i == 2:
		a=pyp.ylabel('Rate/Bin (Hz)', fontsize=18, y=1.1)
		a=pyp.xlabel(r'Radius$^{2}$ (cm$^{2}$)', fontsize=18, x=1.1)
	if i == 2:
		a=pyp.legend(loc = 'upper right', bbox_to_anchor = (1.9, 1.1))

pyp.show()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Load Dev's S1 spectum since we'll need it to apply radial cuts and to generate the BG model PDF

bes1left_detail, bes1right_detail, bins1rate_detail, bins1percenterr_detail = n.loadtxt('/home/tomaszbi/lux/studies/160623_run4coincidenceBack/20160629_s1_rate_fineBins.csv')
bes1mid_detail = (bes1left_detail + bes1right_detail)/2.
bes1width_detail = bes1right_detail - bes1left_detail

bes1leftini, bes1rightini, bins1rateiniini, bins1percenterrini = n.loadtxt('/home/tomaszbi/lux/studies/160623_run4coincidenceBack/20160623_s1_rate.csv')
bes1midini = (bes1leftini+bes1rightini)/2.
bes1widthini = bes1rightini[0]-bes1leftini[0]
bins1rateini = bins1rateiniini
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
"""# correct for the change from spike count to area during the s1 spectrum generation

# S1 histogram of area and spikes
#bes = n.arange(1,10,.1)
bes = n.arange(bes1leftini[0],bes1rightini[-1]+bes1widthini,bes1widthini)
i=2

histdata_area, nothing = n.histogram(rq[i]['pulse_area_phe'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]][fininds], bins=bes)
histdata_area = n.float64(histdata_area) / livetimes[i]

histdata_spikes, nothing = n.histogram(pulse_cspike_phe[i][fininds], bins=bes)
histdata_spikes = n.float64(histdata_spikes) / livetimes[i]

a=pyp.loglog()
a=pyp.plot(bes[1:], histdata_area, 'b-', drawstyle='steps-pre')
a=pyp.plot(bes[:-1], histdata_area, 'b-', drawstyle='steps-post')
a=pyp.plot(bes[1:], histdata_spikes, 'g-', drawstyle='steps-pre')
a=pyp.plot(bes[:-1], histdata_spikes, 'g-', drawstyle='steps-post')
pyp.show()
"""

"""
spike_area_ratio = n.float64(histdata_spikes_temp)/histdata_area_temp
a=pyp.plot(bes[1:], spike_area_ratio, 'g-', drawstyle='steps-pre')
a=pyp.plot(bes[:-1], spike_area_ratio, 'g-', drawstyle='steps-post')
pyp.show()

# apply this ratio to the read in S1 spectrum and correct for the rates
normrate = bins1rateiniini[:10].sum()
bins1rateini = bins1rateiniini * spike_area_ratio
bins1rateini *= normrate/bins1rateini[:10].sum()

a=pyp.loglog()
a=pyp.plot(bes1midini,bins1rateiniini)
a=pyp.plot(bes1midini,bins1rateini)
pyp.show()
"""
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# Remove the hump at 25phe.
nonhumpinds = (bes1leftini <= 11) | ((bes1leftini >= 41) & (bes1leftini <= 46)) | (bes1leftini >= 65)
#bes1left, bes1right = bes1leftini[nonhumpinds], bes1rightini[nonhumpinds]
#bes1mid = bes1midini[nonhumpinds]
"""bes1left, bes1right = bes1leftini[nonhumpinds]+0.5, bes1rightini[nonhumpinds]+0.5
bes1mid = bes1midini[nonhumpinds]+0.5"""
bes1left, bes1right = bes1leftini[nonhumpinds], bes1rightini[nonhumpinds]
bes1mid = bes1midini[nonhumpinds]
bins1rate, bins1percenterr = bins1rateini[nonhumpinds], bins1percenterrini[nonhumpinds]
bes1width = bes1right[0]-bes1left[0]


# interpolate and assign the edge a value.
bins1ratio = bins1rate/bes1width
bins1ratio_detail = bins1rate_detail/bes1width_detail
#bins1rateinterp = interp1d(bes1mid, bins1ratio, bounds_error=False, fill_value=n.nan)

"""# We are going to merge a fine resolution and a course histogram so before interpolating.
cutind = n.searchsorted(bes1right, 3)
cutind_detail = n.searchsorted(bes1right_detail, 3)
tempbeleft = n.append(bes1left_detail[:cutind_detail+1], bes1left[cutind+1:]) 
tempberight = n.append(bes1right_detail[:cutind_detail+1], bes1right[cutind+1:])
tempbem = (tempbeleft + tempberight)/2.
temprate = n.append(bins1rate_detail[:cutind_detail+1] * bins1rate[:cutind+1].sum()/bins1rate_detail[:cutind_detail+1].sum(), bins1rate[cutind+1:])
mergedbewith = tempberight-tempbeleft
bins1rateinterp = interp1d(tempbem, temprate/mergedbewith, bounds_error=False, fill_value=n.nan)"""


merge_spike_count = True

if merge_spike_count:
	#@@@@@@@@@@@ We're going to replace the low area spectrum with one from RQ
	bes = n.arange(bes1leftini[0],bes1rightini[-1]+bes1widthini,bes1widthini)
	histdata_spikes = n.zeros(bes.size - 1, dtype=n.float64)
	for i in range(rqnum):
		fininds = n.where(n.isfinite(pulse_cspike_phe[i]))[0]
		histdata_spikes_temp, nothing = n.histogram(pulse_cspike_phe[i][fininds], bins=bes)
		histdata_spikes += histdata_spikes_temp
	histdata_spikes /= n.sum(livetimes)
	# merge
	cutind = n.searchsorted(bes1right, 11)
	temprate = n.append(histdata_spikes[:cutind+1] * bins1rate[:cutind+1].sum()/histdata_spikes[:cutind+1].sum(), bins1rate[cutind+1:])
	bins1rateinterp = interp1d(bes1mid, temprate/bes1width, bounds_error=False, fill_value=n.nan)
else:
	# use the input S1 area spectrum instead
	bins1rateinterp = interp1d(bes1mid, bins1rate/bes1width, bounds_error=False, fill_value=n.nan)

bes1interp = n.arange(0.5,101, 0.1)
bes1interpmid = (bes1interp[1:] + bes1interp[:-1])/2
bes1interpwidth = bes1interp[1]-bes1interp[0]

# make the interpolated array
finebins1rate = bins1rateinterp(bes1interpmid)
#finebins1rate[bes1interpmid < bes1mid_detail[0]] = 0
finebins1rate[bes1interpmid < bes1mid[0]] = 0
finebins1rate[bes1interpmid > bes1mid[-1]] = bins1ratio[-1]
finebins1rate = finebins1rate * bins1rate.sum()/finebins1rate.sum() * bes1width/bes1interpwidth

############## 
# compare Dev's spectrum to the RQ S1s

bes = n.unique(n.append(bes1leftini, bes1rightini))
i=0

histdata_area, nothing = n.histogram(rq[i]['pulse_area_phe'][s1inds_iso[i]][passed_s1_iso_bad_area_flags[i]], bins=bes)
histdata_area = n.float64(histdata_area) / livetimes[i]

histdata_spikes, nothing = n.histogram(pulse_cspike_phe[i], bins=bes)
histdata_spikes = n.float64(histdata_spikes) / livetimes[i]

a=pyp.loglog()
a=pyp.plot(bes1rightini, bins1rateini/bes1widthini, color='b', lw=1.5, drawstyle='steps-pre', label='DATs')
a=pyp.plot(bes1leftini, bins1rateini/bes1widthini, color='b', lw=1.5, drawstyle='steps-post')

a=pyp.plot(bes1right, bins1rate/bes1width, color='g', ls='--', lw=2, drawstyle='steps-pre', label='DATs, No Hump')
a=pyp.plot(bes1left, bins1rate/bes1width, color='g', ls='--', lw=2, drawstyle='steps-post')

a=pyp.plot(bes[1:], histdata_area/bes1widthini, color='orange', lw=1.5, drawstyle='steps-pre', label='RQ, Area')
a=pyp.plot(bes[:-1], histdata_area/bes1widthini, color='orange', lw=1.5, drawstyle='steps-post')

a=pyp.plot(bes[1:], histdata_spikes/bes1widthini, color='c', lw=1.5, drawstyle='steps-pre', label = 'RQ, Overlap-Corrected Spikes')
a=pyp.plot(bes[:-1], histdata_spikes/bes1widthini, color='c', lw=1.5, drawstyle='steps-post')

a=pyp.plot(bes1interp[1:], finebins1rate, color='r', drawstyle='steps-pre', label='Model Used')
a=pyp.plot(bes1interp[:-1], finebins1rate, color='r', drawstyle='steps-post')
a=pyp.xlabel('S1 (phd)', fontsize=18)
a=pyp.ylabel('Rate/Pulse (Hz/phd)', fontsize=18)
a=pyp.xlim([1,100])
a=pyp.legend()
pyp.show()




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#$$$ Generate an S1 "simulation" from the histogram.
scale_factor_badAreaOnly = 10
scale_factor = 100
#simnums = n.array([scale_factor * s2inds_iso[i][0].size for i in range(rqnum)])
simnums_badAreaOnly = n.array([scale_factor_badAreaOnly * passed_s2_iso_bad_area_flags[i].sum() for i in range(rqnum)])
simnums = n.array([scale_factor * total_s2_flags[i].sum() for i in range(rqnum)])
s1sims_badAreaOnly = []
s1sims = []
n.random.seed(42)
for i in range(rqnum):
	s1be = bes1interp
	s1bm = bes1interpmid
	s1binned_hz = finebins1rate
	s1frac=s1binned_hz/s1binned_hz.sum()
	# Define the continous bounds that we'll use to cut the randomly drawn numbers.
	# Define in such a way that we can use n.searchsorted() to look for which bin each
	# random draw belongs to.
	s1integfrac=n.array([s1frac[:j+1].sum() for j in range(s1frac.size-1)]+[1])
	# Assign simulated counts to each bin.
	####
	s1randbins = n.searchsorted(s1integfrac, n.random.rand(simnums[i]))
	# create a distribution with continous values of S1
	s1bestep = s1be[1]-s1be[0]
	s1sims.append(s1bm[s1randbins] + (n.random.rand(simnums[i])-.5)*s1bestep)
	####
	s1randbins = n.searchsorted(s1integfrac, n.random.rand(simnums_badAreaOnly[i]))
	# create a distribution with continous values of S1
	s1bestep = s1be[1]-s1be[0]
	s1sims_badAreaOnly.append(s1bm[s1randbins] + (n.random.rand(simnums_badAreaOnly[i])-.5)*s1bestep)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# get actual rates in the region where we have data (for S2) and the histogram (for S1)

official_s2_rate_badAreaOnly =  [passed_s2_iso_bad_area_flags[i].sum()/livetimes[i] for i in range(rqnum)]
official_s2_rate =  [total_s2_flags[i].sum()/livetimes[i] for i in range(rqnum)]

official_s1_rate = finebins1rate.sum()*bes1interpwidth

# coincidence rate per day
total_coincident_rate_badAreaOnly = [(official_s1_rate * official_s2_rate_badAreaOnly[i] * (300-40) * 1e-6*24.*3600.) for i in range(rqnum)]

total_coincident_rate = [(official_s1_rate * official_s2_rate[i] * (300-40) * 1e-6*24.*3600.) for i in range(rqnum)]

# compute the factors needed for normalization of the PDF)
multip_factors_badAreaOnly = [simnums_badAreaOnly[i] / total_coincident_rate_badAreaOnly[i] for i in range(rqnum)]
multip_factors = [simnums[i] / total_coincident_rate[i] for i in range(rqnum)]

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# generate drift times
n.random.seed(242)
dtmax = 300.
dtmin=40.
dtrange = dtmax - dtmin

dts = []
dts_badAreaOnly = []
for i in range(rqnum):
	dts.append(n.random.rand(simnums[i]) * dtrange + dtmin)
	dts_badAreaOnly.append(n.random.rand(simnums_badAreaOnly[i]) * dtrange + dtmin)

# Create paired lists
simmed_s2_inds_badAreaOnly = []
simmed_s2_inds = []
for i in range(rqnum):
	# bad Area only
	# rows
	rowinds = s2inds_iso[i][0][passed_s2_iso_bad_area_flags[i]]
	# columns
	colinds = s2inds_iso[i][1][passed_s2_iso_bad_area_flags[i]]
	randomizing = n.tile(n.arange(passed_s2_iso_bad_area_flags[i].sum()),scale_factor_badAreaOnly)
	n.random.shuffle(randomizing)
	simmed_s2_inds_badAreaOnly.append( (rowinds[randomizing], colinds[randomizing]) )
	#### S2 quality cuts
	# rows
	rowinds = s2inds_iso[i][0][total_s2_flags[i]]
	# columns
	colinds = s2inds_iso[i][1][total_s2_flags[i]]
	randomizing = n.tile(n.arange(total_s2_flags[i].sum()),scale_factor)
	n.random.shuffle(randomizing)
	simmed_s2_inds.append( (rowinds[randomizing], colinds[randomizing]) )


##############################################################################
# some quick validation plots
"""
bes = 10**n.linspace(2,4,40)
i=0

a=pyp.loglog()
histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][simmed_s2_inds_badAreaOnly[i]], bins=bes)
histdata = n.float64(histdata) / livetimes[i] / scale_factor
a=pyp.plot(bes[1:], histdata, 'b', lw=2, drawstyle='steps-pre', label='Isolated S2')
a=pyp.plot(bes[:-1], histdata, 'b', lw=2, drawstyle='steps-post')
pyp.show()




i=0
a=pyp.hist(r_cm_sq[i][simmed_s2_inds_badAreaOnly[i]])
a=pyp.show()

a=pyp.plot(r_cm_sq[i][simmed_s2_inds_badAreaOnly[i]], dts_badAreaOnly[i], '.')
pyp.show()
"""
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# apply the radial cut to see how efficient we are.


# Load in the radial fiducial cut info.
bigtimebin, dtlowbe, dthighbe, philowbe, phihighbe, radcutini = n.loadtxt('/home/tomaszbi/lux/studies/160623_run4coincidenceBack/myanalysis/wallposition_run4.txt', delimiter=',', unpack=True)

wall_offset = 3.0	# We don't know at this time if we'll go with 4, 3.5 or 3. Hence I'm going with 3.5cm - this introduces a ~9% bias in either direction if we end up going with 3cm or 4cm. I can live with that.
# On July 1st, we decided to go with a 3cm wall
radcut = radcutini - wall_offset

simfidflags_badAreaOnly = []
simfidflags = []
for i in range(rqnum):
	simfidflags_badAreaOnly.append(n.zeros(simnums_badAreaOnly[i], dtype=n.bool))
	simfidflags.append(n.zeros(simnums[i], dtype=n.bool))
	if i == 0: useinds = n.where(bigtimebin == 0)[0]
	elif i == 1: useinds = n.where(bigtimebin == 2)[0]
	elif i == 2: useinds = n.where(bigtimebin == 3)[0]
	usedtlow, usedthigh = -1, -1
	for j in range(useinds.size):
		if dtlowbe[useinds[j]] != usedtlow:
			usedtlow = dtlowbe[useinds[j]]
			usedthigh = dthighbe[useinds[j]]
			dtflags_badAreaOnly = (dts_badAreaOnly[i] >= usedtlow) & (dts_badAreaOnly[i] < usedthigh)
			dtflags = (dts[i] >= usedtlow) & (dts[i] < usedthigh)
		if dtflags_badAreaOnly.sum() > 0:
			spaceflags = (phi[i][simmed_s2_inds_badAreaOnly[i]][dtflags_badAreaOnly] >= philowbe[useinds[j]]) & (phi[i][simmed_s2_inds_badAreaOnly[i]][dtflags_badAreaOnly] < phihighbe[useinds[j]]) & (r_cm[i][simmed_s2_inds_badAreaOnly[i]][dtflags_badAreaOnly] < radcut[useinds[j]])
			if spaceflags.sum() == 0: continue
			simfidflags_badAreaOnly[i][dtflags_badAreaOnly] = simfidflags_badAreaOnly[i][dtflags_badAreaOnly] | spaceflags
		if dtflags.sum() > 0:
			spaceflags = (phi[i][simmed_s2_inds[i]][dtflags] >= philowbe[useinds[j]]) & (phi[i][simmed_s2_inds[i]][dtflags] < phihighbe[useinds[j]]) & (r_cm[i][simmed_s2_inds[i]][dtflags] < radcut[useinds[j]])
			if spaceflags.sum() == 0: continue
			simfidflags[i][dtflags] = simfidflags[i][dtflags] | spaceflags

for i in range(rqnum):
	print 'Radial cut on %d kept %f events' %(i, float(simfidflags_badAreaOnly[i].sum())/simnums_badAreaOnly[i])

# with 4cm offset from the wall
#Radial cut on 0 kept 0.146178 events
#Radial cut on 1 kept 0.137381 events
#Radial cut on 2 kept 0.102583 events

# with 3.5cm offset from the wall
#Radial cut on 0 kept 0.159011 events
#Radial cut on 1 kept 0.149755 events
#Radial cut on 2 kept 0.111818 events


# with 3cm offset from the wall
#Radial cut on 0 kept 0.171965 events
#Radial cut on 1 kept 0.162501 events
#Radial cut on 2 kept 0.121966 events


for i in range(rqnum):
	print 'Radial cut on %d kept %f events' %(i, float(simfidflags[i].sum())/simnums[i])

# with 4cm offset from the wall
#Radial cut on 0 kept 0.183942 events
#Radial cut on 1 kept 0.170827 events
#Radial cut on 2 kept 0.174235 events

# with 3.5cm offset from the wall
#Radial cut on 0 kept 0.199795 events
#Radial cut on 1 kept 0.186098 events
#Radial cut on 2 kept 0.189592 events

# with 3cm offset from the wall
#Radial cut on 0 kept 0.217379 events
#Radial cut on 1 kept 0.202466 events
#Radial cut on 2 kept 0.205896 events


"""
phibe = n.linspace(-n.pi, n.pi, 29)
dtbe = n.linspace(40, 300, 27)
dtlowbe = n.repeat(dtbe[:-1], 28)
dthighbe = n.repeat(dtbe[1:], 28)
philowbe = n.tile(phibe[:-1], 26)
phihighbe = n.tile(phibe[1:], 26)

simfidflags_badAreaOnly = []
simfidflags = []
for i in range(rqnum):
	simfidflags_badAreaOnly.append(n.zeros(simnums_badAreaOnly[i], dtype=n.bool))
	simfidflags.append(n.zeros(simnums[i], dtype=n.bool))
	if i == 0: radcut = n.loadtxt('fiducialRMax_TB1.txt')
	elif i == 1: radcut = n.loadtxt('fiducialRMax_TB3.txt')
	elif i == 2: radcut = n.loadtxt('fiducialRMax_TB4.txt')
	radcut += 
	usedtlow, usedthigh = -1, -1
	for j in range(radcut.size):
		if dtlowbe[j] != usedtlow:
			usedtlow = dtlowbe[j]
			usedthigh = dthighbe[j]
			dtflags_badAreaOnly = (dts_badAreaOnly[i] >= usedtlow) & (dts_badAreaOnly[i] < usedthigh)
			dtflags = (dts[i] >= usedtlow) & (dts[i] < usedthigh)
		if dtflags_badAreaOnly.sum() > 0: 
			spaceflags = (phi[i][simmed_s2_inds_badAreaOnly[i]][dtflags_badAreaOnly] >= philowbe[j]) & (phi[i][simmed_s2_inds_badAreaOnly[i]][dtflags_badAreaOnly] < phihighbe[j]) & (r_cm[i][simmed_s2_inds_badAreaOnly[i]][dtflags_badAreaOnly] < radcut[j])
			if spaceflags.sum() == 0: continue
			simfidflags_badAreaOnly[i][dtflags_badAreaOnly] = simfidflags_badAreaOnly[i][dtflags_badAreaOnly] | spaceflags
		if dtflags.sum() > 0: 
			spaceflags = (phi[i][simmed_s2_inds[i]][dtflags] >= philowbe[j]) & (phi[i][simmed_s2_inds[i]][dtflags] < phihighbe[j]) & (r_cm[i][simmed_s2_inds[i]][dtflags] < radcut[j])
			if spaceflags.sum() == 0: continue
			simfidflags[i][dtflags] = simfidflags[i][dtflags] | spaceflags


for i in range(rqnum):
	print 'Radial cut on %d kept %f events' %(i, float(simfidflags_badAreaOnly[i].sum())/simnums_badAreaOnly[i])

#Radial cut on 0 kept 0.146178 events
#Radial cut on 1 kept 0.137381 events
#Radial cut on 2 kept 0.102583 events

for i in range(rqnum):
	print 'Radial cut on %d kept %f events' %(i, float(simfidflags[i].sum())/simnums[i])

#Radial cut on 0 kept 0.183942 events
#Radial cut on 1 kept 0.170827 events
#Radial cut on 2 kept 0.174235 events
"""

"""
# plot r-dt plane - bad area only cut
i=0
a=pyp.plot(r_cm_sq[i][simmed_s2_inds_badAreaOnly[i]], dts_badAreaOnly[i], '.')
a=pyp.plot(r_cm_sq[i][simmed_s2_inds_badAreaOnly[i]][simfidflags_badAreaOnly[i]], dts_badAreaOnly[i][simfidflags_badAreaOnly[i]], 'o', mfc='None', mec='r', mew=1.5, ms=7)
pyp.show()

# plot r-dt plane - full cuts
i=0
a=pyp.plot(r_cm_sq[i][simmed_s2_inds[i]], dts[i], '.', label='Simulated Events')
a=pyp.plot(r_cm_sq[i][simmed_s2_inds[i]][simfidflags[i]], dts[i][simfidflags[i]], 'o', mfc='None', mec='r', mew=1.5, ms=7, label='Passed Radial Cut')
a=pyp.xlim([0,21**2])
a=pyp.ylim([35, 305])
a=pyp.legend()
a=pyp.xlabel(r'Radius$^{2}$ (cm$^{2}$)', fontsize=18)
a=pyp.ylabel(r'Drift Time ($\mu$s)', fontsize=18)
pyp.show()


# plot x-y plane - bad area only cut
i=0
a=pyp.plot(rq[i]['x_cm'][simmed_s2_inds_badAreaOnly[i]], rq[i]['y_cm'][simmed_s2_inds_badAreaOnly[i]], '.')
a=pyp.plot(rq[i]['x_cm'][simmed_s2_inds_badAreaOnly[i]][simfidflags_badAreaOnly[i]], rq[i]['y_cm'][simmed_s2_inds_badAreaOnly[i]][simfidflags_badAreaOnly[i]], 'o', mfc='None', mec='r', mew=1.5, ms=7)
pyp.show()

# plot x-y plane - full cuts
i=0
a=pyp.plot(rq[i]['x_cm'][simmed_s2_inds[i]], rq[i]['y_cm'][simmed_s2_inds[i]], '.', label='Simulated Events')
a=pyp.plot(rq[i]['x_cm'][simmed_s2_inds[i]][simfidflags[i]], rq[i]['y_cm'][simmed_s2_inds[i]][simfidflags[i]], 'o', mfc='None', mec='r', mew=1.5, ms=7, label='Passed Radial Cut')
a=pyp.xlim([-25,25])
a=pyp.ylim([-25,25])
a=pyp.legend()
a=pyp.xlabel('X Position (cm)', fontsize=18)
a=pyp.ylabel('Y Position (cm)', fontsize=18)
pyp.show()


# compare r-dt for all 3 samples
i=0
a=pyp.plot(r_cm_sq[i][simmed_s2_inds_badAreaOnly[i]][simfidflags_badAreaOnly[i]], dts_badAreaOnly[i][simfidflags_badAreaOnly[i]], 'b,')
i=1
a=pyp.plot(r_cm_sq[i][simmed_s2_inds_badAreaOnly[i]][simfidflags_badAreaOnly[i]], dts_badAreaOnly[i][simfidflags_badAreaOnly[i]], 'g,')
i=2
a=pyp.plot(r_cm_sq[i][simmed_s2_inds_badAreaOnly[i]][simfidflags_badAreaOnly[i]], dts_badAreaOnly[i][simfidflags_badAreaOnly[i]], 'r,')
pyp.show()
"""

##############
# Compute the drift time corrections to alter the S2 size
dt_corr_factor_slope = (1.3-1.)/(300.-0)
dt_corr_factors_badAreaOnly = [(dts_badAreaOnly[i]*dt_corr_factor_slope+1) for i in range(rqnum)]
dt_corr_factors = [(dts[i]*dt_corr_factor_slope+1) for i in range(rqnum)]

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Compute the final cut around the bands
band_cut_kept_flags_badAreaOnly = []
band_cut_kept_flags = []
final_s2Csim_badAreaOnly = []
final_s2Csim = []
final_s2Csim_nonBandCut = []
final_s1sim_badAreaOnly = []
final_s1sim = []
final_s1sim_nonBandCut = []
for i in range(rqnum):
	temps1 = s1sims[i]
	temps2 = rq[i]['pulse_area_phe'][simmed_s2_inds[i]] * dt_corr_factors[i]
	ltemps1 = n.log10(temps1)
	ltemps2 = n.log10(temps2)
	band_cut_kept_flags.append(
		(ltemps2 < 3.326*(temps1**(-0.084)) + ltemps1) &
		(ltemps2 > 0.8823*(temps1**0.096) + ltemps1)
		)
	final_s2Csim.append(temps2[simfidflags[i] & band_cut_kept_flags[i]])
	final_s1sim.append(temps1[simfidflags[i] & band_cut_kept_flags[i]])
	final_s2Csim_nonBandCut.append(temps2[simfidflags[i] ])
	final_s1sim_nonBandCut.append(temps1[simfidflags[i] ])
	###
	temps1 = s1sims_badAreaOnly[i]
	temps2 = rq[i]['pulse_area_phe'][simmed_s2_inds_badAreaOnly[i]] * dt_corr_factors_badAreaOnly[i]
	ltemps1 = n.log10(temps1)
	ltemps2 = n.log10(temps2)
	band_cut_kept_flags_badAreaOnly.append(
		(ltemps2 < 3.326*(temps1**(-0.084)) + ltemps1) &
		(ltemps2 > 0.8823*(temps1**0.096) + ltemps1)
		)
	final_s2Csim_badAreaOnly.append(temps2[simfidflags_badAreaOnly[i] & band_cut_kept_flags_badAreaOnly[i]])
	final_s1sim_badAreaOnly.append(temps1[simfidflags_badAreaOnly[i] & band_cut_kept_flags_badAreaOnly[i]])

#Upper 3-Sigma-Above-ER-Band Cut: log10(S2zc) < 3.326*(S1^(-0.084)) + log10(S1)
#Lower 7-Sigma-Below-NR-Band Cut: log10(S2zc) > 0.8823*(S1^0.096) + log10(S1)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




# Get the spectrum of final simulated S2s along with the various cuts




bes = 10**n.linspace(2,4,40)

a=pyp.figure(1, figsize=(10,8))
for i in range(rqnum):
	a=pyp.suptitle('Isolated S2 Distributions', fontsize=18)
	a=pyp.subplot(2,2,i+1)
	a=pyp.grid(which='both')
	a=pyp.loglog()
	a=pyp.title(dsets[i])
	
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][simmed_s2_inds[i]], bins=bes)
	#histdata, nothing = n.histogram(final_s2Csim[i], bins=bes)
	histdata = n.float64(histdata) / scale_factor/livetimes[i] #multip_factors[i]
	a=pyp.plot(bes[1:], histdata, color='r', lw=5, drawstyle='steps-pre', label='S2 Selection Cuts')
	a=pyp.plot(bes[:-1], histdata, color='r', lw=5, drawstyle='steps-post')
	
	histdata, nothing = n.histogram(final_s2Csim[i], bins=bes)
	histdata = n.float64(histdata) / scale_factor/livetimes[i] #multip_factors[i]
	a=pyp.plot(bes[1:], histdata, color='r', lw=5, drawstyle='steps-pre', label='S2 + Band + Fid Cuts')
	a=pyp.plot(bes[:-1], histdata, color='r', lw=5, drawstyle='steps-post')
	
	
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, 'b', lw=2, drawstyle='steps-pre', label='All')
	a=pyp.plot(bes[:-1], histdata, 'b', lw=2, drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, 'g', lw=2, drawstyle='steps-pre', label='Bad Area Cut')
	a=pyp.plot(bes[:-1], histdata, 'g', lw=2, drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & sigThr_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, '--', color='c', drawstyle='steps-pre', label='Bad Area Cut + Sigma Thrsh.')
	a=pyp.plot(bes[:-1], histdata, '--', color='c',  drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & aftSigRatio_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, '--', color='orange', drawstyle='steps-pre', label='Bad Area Cut + AFT-Sigma Ratio')
	a=pyp.plot(bes[:-1], histdata, '--', color='orange',  drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & chi2cut_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, '--', color='m', drawstyle='steps-pre', label=r'Bad Area Cut + $\chi^{2}$')
	a=pyp.plot(bes[:-1], histdata, '--', color='m',  drawstyle='steps-post')
	#
	histdata, nothing = n.histogram(rq[i]['pulse_area_phe'][s2inds_iso[i]][passed_s2_iso_bad_area_flags[i] & sigThr_flags[i] & aftSigRatio_flags[i] & chi2cut_flags[i]], bins=bes)
	histdata = n.float64(histdata) / livetimes[i]
	a=pyp.plot(bes[1:], histdata, '-', color='k', lw=1.5, drawstyle='steps-pre', label=r'Bad Area Cut + All S2 Cuts')
	a=pyp.plot(bes[:-1], histdata, '-', color='k', lw=1.5,  drawstyle='steps-post')
	#
	if i == 2:
		a=pyp.ylabel('Rate/Bin (Hz)', fontsize=18, y=1.1)
		a=pyp.xlabel('Raw S2 Area (phd)', fontsize=18, x=1.1)
	if i == 2:
		a=pyp.legend(loc = 'upper right', bbox_to_anchor = (1.9, 1.1))

pyp.show()

############# 
# write out the data so that we can use is somewhere else

# get the total livetime
tl = (n.array(livetimes)).sum()*scale_factor


# ------ The S2 areas after all data quality cuts (includes drift time correction)
# store the livetime
newtxt = "%.9e\n" %(tl)
# go through the S2s
cc=0
for i in range(rqnum):
    temps2 = rq[i]['pulse_area_phe'][simmed_s2_inds[i]] * dt_corr_factors[i]
    for j in range(temps2.size):
        newtxt += '%.4e\n' %(temps2[j])
        cc += 1

f=open('isolatedS2s_livetime_areas_allS2QualityCuts_180513.txt', 'w')
f.write(newtxt)
f.close()


# ------ The S2 areas after all data quality and fiducial cuts. band cut used
# get the total livetime
tl = (n.array(livetimes)).sum()*scale_factor
# store the livetime
newtxt = "%.9e\n" %(tl)
# go through the S2s
cc=0
for i in range(rqnum):
    for j in range(final_s2Csim[i].size):
        newtxt += '%.4e\n' %(final_s2Csim[i][j])
        cc+=1

f=open('isolatedS2s_livetime_areas_allS2QualityCuts-fidCuts_180517.txt', 'w')
f.write(newtxt)
f.close()


# ------ The S2 areas after all data quality and fiducial cuts. non-band cut used
# get the total livetime
tl = (n.array(livetimes)).sum()*scale_factor
# store the livetime
newtxt = "%.9e\n" %(tl)
# go through the S2s
cc=0
for i in range(rqnum):
    for j in range(final_s2Csim_nonBandCut[i].size):
        newtxt += '%.4e\n' %(final_s2Csim_nonBandCut[i][j])
        cc+=1

f=open('isolatedS2s_livetime_areas_allS2QualityCuts-fidCuts-nonBandCut_180517.txt', 'w')
f.write(newtxt)
f.close()



# ------ The S1 areas after all data quality and fiducial cuts. band-cut used
# get the total livetime
tl = (n.array(livetimes)).sum()*scale_factor
# store the livetime
newtxt = "%.9e\n" %(tl)
# go through the S1s
cc=0
for i in range(rqnum):
    for j in range(final_s1sim[i].size):
        newtxt += '%.4e\n' %(final_s1sim[i][j])
        cc+=1

f=open('isolatedS1s_livetime_areas_allS2QualityCuts-fidCuts_180517.txt', 'w')
f.write(newtxt)
f.close()


tl = (n.array(livetimes)).sum()*scale_factor
# store the livetime
newtxt = "%.9e\n" %(tl)
# go through the S1s
cc=0
for i in range(rqnum):
    for j in range(final_s1sim_nonBandCut[i].size):
        newtxt += '%.4e\n' %(final_s1sim_nonBandCut[i][j])
        cc+=1

f=open('isolatedS1s_livetime_areas_allS2QualityCuts-fidCuts-nonBandCut_180517.txt', 'w')
f.write(newtxt)
f.close()

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






# Get the PDF


class my_format(Formatter):
	"""
	"""
	def __init__(self, extrapobj, log=False):
		self.extrapobj = extrapobj
		self.log = log
		Formatter.__init__(self)
	def __call__(self, x, pos=None):
		if self.log:	# log
			newticklabel = '%.1e' %(10**self.extrapobj(x))
		else: # linear
			newticklabel = '%.1f' %(self.extrapobj(x))
		return newticklabel



# Get the approximate NR band
nrs1 = n.array([2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50])
nrmu = n.array([2.3679,2.2185,2.1569,2.1159,2.0748,2.0493,2.0119,1.9986,\
                1.9697,1.9439,1.9236,1.9057,1.9023,1.8642,1.8685,1.8657,\
                1.8506,1.8463,1.8322,1.8215,1.8181,1.8068,1.7833,1.7555,1.7559])
nrlo = nrmu-1.28*n.array([0.216,0.1935,0.1611,0.134,0.1281,0.1214,0.1136,0.1101,0.0993,\
                            0.1151,0.1163,0.0954,0.1052,0.1038,0.0939,0.0785,0.0827,\
                            0.0809,0.0803,0.0884,0.0837,0.08,0.0809,0.0686,0.0755])
nrhi = nrmu+1.28*n.array([0.216,0.1935,0.1611,0.134,0.1281,0.1214,0.1136,0.1101,0.0993,\
                            0.1151,0.1163,0.0954,0.1052,0.1038,0.0939,0.0785,0.0827,\
                            0.0809,0.0803,0.0884,0.0837,0.08,0.0809,0.0686,0.0755])

# extrapolate the band after getting it into S1-log10(S2) space
nlog10S1 = n.log10(nrs1)
ls2bandextrap10_nr = extrap1d.extrap1d(nrs1, nrlo+nlog10S1 )
ls2bandextrap_nr = extrap1d.extrap1d(nrs1, nrmu+nlog10S1 )
ls2bandextrap90_nr = extrap1d.extrap1d(nrs1, nrhi+nlog10S1 )



# bin the data according to the analysis binning
s1bes = n.linspace(1,50,50)
s2bes = n.linspace(2,4,41)
#
#s1bes = n.linspace(1,50,10)
#s2bes = n.linspace(2,4,10)
s1bm = (s1bes[1:]+s1bes[:-1])/2.

# prepare a matrix for normalizing the PDFs
pdfnormmattemp = n.meshgrid(s2bes[1:]-s2bes[:-1], s1bes[1:]-s1bes[:-1])	


xmat_badAreaOnly, ymat_badAreaOnly, counts_badAreaOnly = [], [], []
lcounts_badAreaOnly = []
pdf_badAreaOnly, lpdf_badAreaOnly = [], []
xmat, ymat, counts = [], [], []
lcounts = []
pdf, lpdf = [], []
for i in range(rqnum):
	# Take the randomized S2 areas that survivied the bad area cut and a fiducial radial cut and apply smearing due to the drift time correction.
	s2sim_badAreaOnly = final_s2Csim_badAreaOnly[i]
	s2sim = final_s2Csim[i]
	s1sim_badAreaOnly = final_s1sim_badAreaOnly[i]
	s1sim = final_s1sim[i]
	#
	# Bin in S1 and S2 space
	temp = hist2d.hist2d(s1bes, s1sim_badAreaOnly, s2bes, n.log10(s2sim_badAreaOnly))
	xmat_badAreaOnly.append(temp[0])
	ymat_badAreaOnly.append(temp[1])
	counts_badAreaOnly.append(temp[2])
	#
	temp = hist2d.hist2d(s1bes, s1sim, s2bes, n.log10(s2sim))
	xmat.append(temp[0])
	ymat.append(temp[1])
	counts.append(temp[2])
	#
	# normalize to counts per day and compute log10 of that
	counts_badAreaOnly[i] = counts_badAreaOnly[i] / multip_factors_badAreaOnly[i]
	lcounts_badAreaOnly.append(n.log10(counts_badAreaOnly[i]))
	counts[i] = counts[i] / multip_factors[i]
	lcounts.append(n.log10(counts[i]))
	#
	# normalize the pdf
	# normailze to 1 to get the actual pdf (instead of events per day)
	pdf_badAreaOnly.append(counts_badAreaOnly[i]/(pdfnormmattemp[0]*pdfnormmattemp[1])/counts_badAreaOnly[i].sum())
	lpdf_badAreaOnly.append(n.log10(pdf_badAreaOnly[i]))
	pdf.append(counts[i]/(pdfnormmattemp[0]*pdfnormmattemp[1])/counts[i].sum())
	lpdf.append(n.log10(pdf[i]))


i=0
# print minima and maxima of various matrices to enable me to set plot ranges
print counts_badAreaOnly[i][n.isfinite(counts_badAreaOnly[i])].min()*375
print counts_badAreaOnly[i][n.isfinite(counts_badAreaOnly[i])].max()*375
print lcounts_badAreaOnly[i][n.isfinite(lcounts_badAreaOnly[i])].min()
print lcounts_badAreaOnly[i][n.isfinite(lcounts_badAreaOnly[i])].max()
print lcounts_badAreaOnly[i][n.isfinite(lcounts_badAreaOnly[i])].min() + n.log10(375)
print lcounts_badAreaOnly[i][n.isfinite(lcounts_badAreaOnly[i])].max() + n.log10(375)

i=0
# print minima and maxima of various matrices to enable me to set plot ranges
print counts[i][n.isfinite(counts[i])].min()*375
print counts[i][n.isfinite(counts[i])].max()*375
print lcounts[i][n.isfinite(lcounts[i])].min()
print lcounts[i][n.isfinite(lcounts[i])].max()
print lcounts[i][n.isfinite(lcounts[i])].min() + n.log10(375)
print lcounts[i][n.isfinite(lcounts[i])].max() + n.log10(375)


avgcounts = n.mean(counts, axis=0)
# print minima and maxima of various matrices to enable me to set plot ranges
print avgcounts[n.isfinite(avgcounts)].min()*375
print avgcounts[n.isfinite(avgcounts)].max()*375
lavgcounts = n.log10(avgcounts)
print lavgcounts[n.isfinite(lavgcounts)].min()+n.log10(375)
print lavgcounts[n.isfinite(lavgcounts)].max()+n.log10(375)




#vmin, vmax = -5, -1.5
#vmin, vmax = -2, 0.7
#vmin, vmax = 0, 4.8
vmin, vmax = 0, .032

normobj = mlib.colors.Normalize(vmin=vmin, vmax=vmax)
#csobs = mlib.cm.ScalarMappable(cmap=mlib.cm.Greys,norm=normobj)
csobs = mlib.cm.ScalarMappable(cmap=mlib.cm.jet,norm=normobj)

xextrap = extrap1d.extrap1d(n.arange(s1bes.size)-0.5, s1bes)
#yextrap = extrap1d.extrap1d(n.arange(ratiobes.size)-0.5, ratiobes)
yextrap = extrap1d.extrap1d(n.arange(s2bes.size)-0.5, s2bes)
# define the inverse extrapolation such that we can overplot bands
ixextrap = extrap1d.extrap1d(s1bes, n.arange(s1bes.size)-0.5)
#iyextrap = extrap1d.extrap1d(ratiobes, n.arange(ratiobes.size)-0.5)
iyextrap = extrap1d.extrap1d(s2bes, n.arange(s2bes.size)-0.5)

a=pyp.figure(1, figsize=(10,8))
a=pyp.subplots_adjust(bottom=0.1,top=0.92, right=0.84, left=0.09)
ax = pyp.subplot(1,1,1)
#a=pyp.imshow(lavgcounts.transpose()+n.log10(375), origin='bottom', interpolation='nearest', cmap=csobs.get_cmap(), norm=normobj, aspect='auto')
a=pyp.imshow(avgcounts.transpose()*375., origin='bottom', interpolation='nearest', cmap=csobs.get_cmap(), norm=normobj, aspect='auto')
cpltx = ixextrap(s1bm)
a=pyp.plot(cpltx,iyextrap(ls2bandextrap_nr(s1bm)),'r-', lw=1.5)#, label='NR Mean')
a=pyp.plot(cpltx,iyextrap(ls2bandextrap10_nr(s1bm)),'r--', lw=1.5)
a=pyp.plot(cpltx,iyextrap(ls2bandextrap90_nr(s1bm)),'r--', lw=1.5)
a=pyp.xlabel(r'$\mathrm{S1\ (phd)}$', fontsize=20)
a=pyp.ylabel(r'$\mathrm{log_{10}[S2]}$', fontsize=20)

# set limits
a=pyp.ylim([iyextrap(2),iyextrap(4)])
a=pyp.xlim([ixextrap(1),ixextrap(50)])

# set ticks
a=ax.xaxis.set_ticks(ixextrap(n.append(1,n.arange(5,50+5,5))))
a=ax.yaxis.set_ticks(iyextrap(n.arange(2,4+.5,.5)))

xformat = my_format(xextrap)
a=ax.xaxis.set_major_formatter(xformat)
yformat = my_format(yextrap)
a=ax.yaxis.set_major_formatter(yformat)

axc = pyp.axes((.85,.1,.015,.82))
cb=mpl.colorbar.ColorbarBase(axc, cmap=csobs.get_cmap(),norm=normobj, format='%0.1f')
#a=axc.yaxis.set_ticks()
if vmax-vmin > 1:
	ticklabels = n.arange(vmin,vmax+.5,.5)
elif vmax-vmin > .35:
	ticklabels = n.arange(vmin,vmax+.05,.05)
elif vmax-vmin > .035:
	ticklabels = n.arange(vmin,vmax+.005,.005)
else:# vmax-vmin > .01:
	ticklabels = n.arange(vmin,vmax+.002,.002)

ticklocs = [normobj(val) for val in ticklabels]
a=axc.yaxis.set_ticks(ticklocs)
a=axc.yaxis.set_ticklabels(ticklabels)
#a=pyp.ylabel(r'$\mathrm{log_{10}[PDF]}$',fontsize=20)
#a=pyp.ylabel(r'$\mathrm{log_{10}[Counts/Bin\ in\ 375\ Days]}$',fontsize=20)
a=pyp.ylabel(r'Counts/Bin in 375 Days',fontsize=20)
for label in ax.get_xticklabels()+ax.get_yticklabels()+axc.get_yticklabels(): label.set_size(18)

pyp.show()



#$$$$$$$$$$$$$$$$$$$$$$

# Print various metrics

rate_region_1 = 0.
rate_region_2 = 0.
rate_160 = 0.
rate_200 = 0.
rate_200raw = 0.
for i in range(rqnum):
	# With bad area cut alone
	s2sim = final_s2Csim[i]
	ls2sim = n.log10(s2sim)
	# number of events within the region of interest
	inregionflags = (final_s1sim[i] > 1) & (final_s1sim[i] < 50) & (ls2sim > 2) & (ls2sim < 4)
	# number of events within the region of ineterest and with S2 larger than 160
	inregionflagsover160 = inregionflags & (s2sim > 160)
	# number of events within the region of ineterest and with S2 larger than 200
	inregionflagsover200 = inregionflags & (s2sim > 200)
	# number of events within the region of ineterest and with S2 larger than 200 before livetime correction
	inregionflagsover200raw = inregionflags & (rq[i]['pulse_area_phe'][simmed_s2_inds[i]][simfidflags[i] & band_cut_kept_flags[i]] > 200)
	rate_region_1 += inregionflags.sum()/multip_factors[i]
	rate_region_2 += counts[i].sum()
	rate_160 += inregionflagsover160.sum()/multip_factors[i]
	rate_200 += inregionflagsover200.sum()/multip_factors[i]
	rate_200raw += inregionflagsover200raw.sum()/multip_factors[i]

rate_region_1 /= rqnum
rate_region_2 /= rqnum
rate_160 /= rqnum
rate_200 /= rqnum
rate_200raw /= rqnum

# total counts from binning
print 'Counts/day in region of interest from binning', rate_region_2

print 'Counts/day within region 1<S1<50, 100<S2<10000', rate_region_1

print 'Counts/day within region with S2>160', rate_160

print 'Counts/day within region with S2>200', rate_200

print 'Counts/day within region with S2>200raw', rate_200raw

# without the cuts around electron and nuclear recoil bands
#Counts/day in region of interest from binning 0.00641382700021
#Counts/day within region 1<S1<50, 100<S2<10000 0.00641382700021
#Counts/day within region with S2>160 0.00590378645159
#Counts/day within region with S2>200 0.00563564260886
#Counts/day within region with S2>200raw 0.00551651261525

# with the cuts around the electron and nuclear recoil bands
#Counts/day in region of interest from binning 0.00452550029683
#Counts/day within region 1<S1<50, 100<S2<10000 0.00452550029683
#Counts/day within region with S2>160 0.0041076047968
#Counts/day within region with S2>200 0.00388041147044
#Counts/day within region with S2>200raw 0.00377875297747


# with the cuts around the electron and nuclear recoil bands and using areas instead of spikes
#Counts/day in region of interest from binning 0.00406454798191
#Counts/day within region 1<S1<50, 100<S2<10000 0.00406454798191
#Counts/day within region with S2>160 0.00371389449919
#Counts/day within region with S2>200 0.00352089566626
#Counts/day within region with S2>200raw 0.00343415543929


#$$$$$$$$$$$$$$$$$$$$$$
# Print out the PDF

avgpdf = n.mean(pdf, axis=0)

"""txt = "# s1 bin edges: "
for be in s1bes: txt+='%.1f ' %(be)

txt += '\n'

txt += "# log10(s2) bin edges: "
for be in s2bes: txt+='%.2f ' %(be)

txt += '\n'
txt += '# S1 bins in columns, log10(S2) bins in rows\n'"""
txt = ''
for i in range(s2bes.size-1):
	for j in range(s1bes.size-1): 
		txt += '%.3e\n' %(avgpdf[j,i])


f=open('run4_randCoincidenceBack_singleLine_final-S1Areas_160714.txt', 'w')
f.write(txt)
f.close()


#$$$$$$$$$$$$$$$$$$$$$$

belowmean = 0.
below90 = 0.
below10 = 0.
for i in range(rqnum):
	# With bad area cut alone
	s2sim_badAreaOnly = final_s2Csim_badAreaOnly[i]
	ls2sim_badAreaOnly = n.log10(s2sim_badAreaOnly)
	# number of events within the region of interest
	inregionflags = (final_s1sim_badAreaOnly[i] > 1) & (final_s1sim_badAreaOnly[i] < 50) & (ls2sim_badAreaOnly > 2) & (ls2sim_badAreaOnly < 4)
	# below band counts
	belowbands2flags_mean_badAreaOnly = inregionflags & (ls2sim_badAreaOnly < ls2bandextrap_nr(final_s1sim_badAreaOnly[i]))
	belowbands2flags_10_badAreaOnly = inregionflags & (ls2sim_badAreaOnly < ls2bandextrap10_nr(final_s1sim_badAreaOnly[i]))
	belowbands2flags_90_badAreaOnly = inregionflags & (ls2sim_badAreaOnly < ls2bandextrap90_nr(final_s1sim_badAreaOnly[i]))
	belowmean += belowbands2flags_mean_badAreaOnly.sum()
	below90 += belowbands2flags_90_badAreaOnly.sum()
	below10 += belowbands2flags_10_badAreaOnly.sum()

print below90/multip_factors_badAreaOnly[i]*375./rqnum

print belowmean/multip_factors_badAreaOnly[i]*375./rqnum

print below10/multip_factors_badAreaOnly[i]*375./rqnum





##############
# With all of the S2 quality cuts


belowmean = 0.
below90 = 0.
below10 = 0.
above90 = 0.
for i in range(rqnum):
	# With bad area cut alone
	s2sim = final_s2Csim[i]
	ls2sim = n.log10(s2sim)
	# number of events within the region of interest
	inregionflags = (final_s1sim[i] > 1) & (final_s1sim[i] < 50) & (ls2sim > 2) & (ls2sim < 4)
	# below band counts
	belowbands2flags_mean = inregionflags & (ls2sim < ls2bandextrap_nr(final_s1sim[i]))
	belowbands2flags_10 = inregionflags & (ls2sim < ls2bandextrap10_nr(final_s1sim[i]))
	belowbands2flags_90 = inregionflags & (ls2sim < ls2bandextrap90_nr(final_s1sim[i]))
	abovebands2flags_90 = inregionflags & (ls2sim > ls2bandextrap90_nr(final_s1sim[i]))
	belowmean += belowbands2flags_mean.sum()
	below90 += belowbands2flags_90.sum()
	below10 += belowbands2flags_10.sum()
	above90 = abovebands2flags_90.sum()



print below90/multip_factors[i]*375./rqnum
print belowmean/multip_factors[i]*375./rqnum
print below10/multip_factors[i]*375./rqnum

# above the band
print above90/multip_factors[i]*375./rqnum


# with the cuts around the electron and nuclear recoil bands
#print below90/multip_factors[i]*375./rqnum
#0.932677132379
#print belowmean/multip_factors[i]*375./rqnum
#0.67883707348
#print below10/multip_factors[i]*375./rqnum
#0.457047712171

## above the band
#print above90/multip_factors[i]*375./rqnum
#0.193678572228

# with the cuts around the electron and nuclear recoil bands and using areas instead of spikes
#>>> print below90/multip_factors[i]*375./rqnum
#0.82439363596
#>>> print belowmean/multip_factors[i]*375./rqnum
#0.601014339636
#>>> print below10/multip_factors[i]*375./rqnum
#0.398577391426
#>>> 
#>>> # above the band
#print above90/multip_factors[i]*375./rqnum
#0.180511310274





##############
# With all of the S2 quality cuts AND the RAW S2 size


belowmean = 0.
below90 = 0.
below10 = 0.
above90 = 0.
for i in range(rqnum):
	# With bad area cut alone
	s2sim = final_s2Csim[i]
	ls2sim = n.log10(s2sim)
	# number of events within the region of interest
	inregionflags = (final_s1sim[i] > 1) & (final_s1sim[i] < 50) & (ls2sim > 2) & (ls2sim < 4) & (rq[i]['pulse_area_phe'][simmed_s2_inds[i]][simfidflags[i] & band_cut_kept_flags[i]] > 200)
	# below band counts
	belowbands2flags_mean = inregionflags & (ls2sim < ls2bandextrap_nr(final_s1sim[i]))
	belowbands2flags_10 = inregionflags & (ls2sim < ls2bandextrap10_nr(final_s1sim[i]))
	belowbands2flags_90 = inregionflags & (ls2sim < ls2bandextrap90_nr(final_s1sim[i]))
	abovebands2flags_90 = inregionflags & (ls2sim > ls2bandextrap90_nr(final_s1sim[i]))
	belowmean += belowbands2flags_mean.sum()
	below90 += belowbands2flags_90.sum()
	below10 += belowbands2flags_10.sum()
	above90 = abovebands2flags_90.sum()



print below90/multip_factors[i]*375./rqnum
print belowmean/multip_factors[i]*375./rqnum
print below10/multip_factors[i]*375./rqnum

# above the band
print above90/multip_factors[i]*375./rqnum



#>>> print below90/multip_factors[i]*375./rqnum
#0.606879577164
#>>> print belowmean/multip_factors[i]*375./rqnum
#0.383500280839
#>>> print below10/multip_factors[i]*375./rqnum
#0.182771151792
#>>> 
#>>> # above the band
#... print above90/multip_factors[i]*375./rqnum
#0.180511310274

print belowmean/multip_factors[i]/rqnum

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# try to estimate errors on the normalization

temp = ((pdf[0]-avgpdf)/n.abs(avgpdf)).flatten()
a=pyp.hist(temp[n.isfinite(temp)])
pyp.show()

temp = ((pdf[1]-avgpdf)/n.abs(avgpdf)).flatten()
a=pyp.hist(temp[n.isfinite(temp)])
pyp.show()

temp = ((pdf[2]-avgpdf)/n.abs(avgpdf)).flatten()
a=pyp.hist(temp[n.isfinite(temp)])
pyp.show()



a=pyp.figure(1)
a=pyp.imshow(pdf[0].transpose(), interpolation='nearest', origin='bottom'); a=pyp.colorbar()
a=pyp.figure(2)
a=pyp.imshow(pdf[1].transpose(), interpolation='nearest', origin='bottom'); a=pyp.colorbar()
a=pyp.figure(3)
a=pyp.imshow(pdf[2].transpose(), interpolation='nearest', origin='bottom'); a=pyp.colorbar()
pyp.show()


