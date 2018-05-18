
import read_root as rr
import numpy as n
import os, sys
import re
#import fitgausshist as fg
import cPickle as cp
import time
#import get_html_data as gethtml
import gc
import glob

#from scipy.optimize import leastsq
#from scipy import stats

#import matplotlib.pyplot as pyp
#import matplotlib as mlib		# for color scale
#from matplotlib import mpl		# for colorbar


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Function definitions
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# RQ settings
branches_1d = ['full_evt_area_phe', 'luxstamp_samples']

branches_2d = ['pulse_classification', 'pulse_start_samples', 'pulse_end_samples', 'pulse_area_phe', 'z_drift_samples', 'aft_t1_samples', 'aft_t0_samples', 'gaus_fit_sigma_samples',
#'selected_s1_s2' - not present in some sets
# mercuty positions
'x_cm', 'y_cm', 'chi2']

branches_3d = ['spike_count', 'peak_area_phe']

# livetime determination
branches_livetime = ['livetime_end_samples', 'livetime_latch_samples']


#rqdirstrs = [['lux10_20141111T0308_cp23714', 'lux10_20141111T0727_cp23713', 'lux10_20141111T1014_cp23712', 'lux10_20141111T1159_cp23442'],
#['lux10_20150701T0404_cp20979', 'lux10_20150701T1351_cp22232', 'lux10_20150701T2215_cp23030', 'lux10_20150702T0610_cp20982'],
#['lux10_20160415T1645_cp21971', 'lux10_20160416T0022_cp21972', 'lux10_20160416T0800_cp21973', 'lux10_20160416T2351_cp21974']]


rqdirstrs = [
	['lux10_20141111T0308_cp23714', 'lux10_20141111T0727_cp23713', 'lux10_20141111T1014_cp23712',
		'lux10_20141111T1159_cp23442', 'lux10_20141111T1419_cp23716', 'lux10_20141111T2226_cp23718',
		'lux10_20141112T0308_cp23722', 'lux10_20141112T0753_cp23724', 'lux10_20141112T1156_cp23728',
		'lux10_20141112T1548_cp22317', 'lux10_20141112T2044_cp23698', 'lux10_20141113T0140_cp23704',
		'lux10_20141113T0638_cp22318', 'lux10_20141113T1046_cp23701', 'lux10_20141113T1453_cp23707',
		'lux10_20141113T1905_cp22319', 'lux10_20141113T2332_cp23710', 'lux10_20141114T0357_cp23711',
		'lux10_20141114T0817_cp23443', 'lux10_20141114T1153_cp23446', 'lux10_20141114T1525_cp23715',
		'lux10_20141114T1934_cp22320', 'lux10_20141115T0001_cp23447', 'lux10_20141115T0433_cp23444',
		'lux10_20141115T0911_cp23719', 'lux10_20141115T1349_cp22321', 'lux10_20141115T1826_cp23445',
		'lux10_20141115T2300_cp23721', 'lux10_20141116T0051_cp22268', 'lux10_20141116T0526_cp23723',
		'lux10_20141116T1000_cp23725', 'lux10_20141116T1435_cp23726', 'lux10_20141116T1916_cp23730',
		'lux10_20141116T2347_cp23448', 'lux10_20141117T0423_cp23731', 'lux10_20141117T0846_cp23732',
		'lux10_20141117T1119_cp23734', 'lux10_20141117T1347_cp23452', 'lux10_20141117T1537_cp23733',
		'lux10_20141117T1932_cp23736'],
	['lux10_20150701T0404_cp20979', 'lux10_20150701T1351_cp22232', 'lux10_20150701T2215_cp23030',
		'lux10_20150702T0610_cp20982', 'lux10_20150702T1405_cp23031', 'lux10_20150702T2152_cp20985',
		'lux10_20150703T0539_cp23032', 'lux10_20150703T1322_cp20987', 'lux10_20150703T2109_cp23033',
		'lux10_20150704T0457_cp20989', 'lux10_20150704T1709_cp23034', 'lux10_20150705T0108_cp20990',
		'lux10_20150705T0151_cp23035', 'lux10_20150705T0959_cp20992', 'lux10_20150705T1805_cp23036',
		'lux10_20150706T0208_cp20995', 'lux10_20150706T0959_cp23037', 'lux10_20150706T1806_cp20997',
		'lux10_20150707T0207_cp23038', 'lux10_20150707T1008_cp20998', 'lux10_20150707T1423_cp23039',
		'lux10_20150707T2220_cp21001', 'lux10_20150708T0615_cp23040', 'lux10_20150708T2037_cp21004',
		'lux10_20150709T0428_cp23042', 'lux10_20150709T1229_cp21005', 'lux10_20150709T1511_cp23043',
		'lux10_20150709T2312_cp21007', 'lux10_20150710T0712_cp23044', 'lux10_20150710T1509_cp21009',
		'lux10_20150710T2313_cp23045', 'lux10_20150711T0712_cp21013', 'lux10_20150711T1516_cp23046',
		'lux10_20150711T2318_cp22226', 'lux10_20150712T0151_cp23047', 'lux10_20150712T0950_cp22229',
		'lux10_20150712T1846_cp23048', 'lux10_20150713T0242_cp21014', 'lux10_20150713T1444_cp21016',
		'lux10_20150713T2242_cp23049'], 
	['lux10_20160415T1645_cp21971', 'lux10_20160416T0022_cp21972', 'lux10_20160416T0800_cp21973',
		'lux10_20160416T2351_cp21974', 'lux10_20160417T0151_cp21975', 'lux10_20160417T0926_cp21976',
		'lux10_20160417T1659_cp21977', 'lux10_20160418T0037_cp21978', 'lux10_20160418T0808_cp21979',
		'lux10_20160418T1346_cp21980', 'lux10_20160418T1718_cp21988', 'lux10_20160419T0052_cp21990',
		'lux10_20160419T0907_cp21991', 'lux10_20160419T0942_cp22000', 'lux10_20160419T1702_cp22005',
		'lux10_20160419T1818_cp22012', 'lux10_20160420T0156_cp22022', 'lux10_20160420T0907_cp22027',
		'lux10_20160420T1642_cp22033', 'lux10_20160421T0020_cp22044', 'lux10_20160421T0750_cp22045',
		'lux10_20160421T1515_cp22046', 'lux10_20160421T2247_cp21989', 'lux10_20160422T0623_cp21999',
		'lux10_20160422T1346_cp22010', 'lux10_20160422T1930_cp22018', 'lux10_20160423T0305_cp22024',
		'lux10_20160423T1038_cp22031', 'lux10_20160423T1805_cp22034', 'lux10_20160424T0142_cp22048',
		'lux10_20160424T0151_cp22061', 'lux10_20160424T0927_cp22062', 'lux10_20160424T1702_cp22063',
		'lux10_20160425T0035_cp22109', 'lux10_20160425T0806_cp22114', 'lux10_20160425T1541_cp22113',
		'lux10_20160425T2310_cp22119', 'lux10_20160426T0637_cp22118', 'lux10_20160426T1312_cp22121',
		'lux10_20160426T2036_cp22122']]

max_files_in_dirs = [None,None,None]
#max_files_in_dirs = [50,50,50]
dsets = ['lux10_201411', 'lux10_201507', 'lux10_201604']

rqdir = '/data/tomaszbi/lux_data/run4/rq'

rqnum = len(rqdirstrs)
# get the rq.root file structure

##############
# Some Tritium and DD events may have snuck in. Cut them out by cutting on luxstamp
livetime_exclude1 = (13110000000000000, 13130000000000000)
livetime_exclude2 = (13301476800000000, 13305000000000000)
livetime_exclude3 = (13336000000000000, 13425000000000000)
livetime_exclude4 = (1.6854e16, n.inf)

# get the bottom indecies
botinds = n.append(n.arange(60,120,1), 121)
topinds = n.append(n.arange(0,60,1), 120)
## remove the bad channels
#botinds = n.delete(botinds, 92-60)
#topinds = n.delete(topinds, [4,31])

# create a timestamp for the run
loctime = time.localtime()
timestamp = '%04d%02d%02dT%02d%02d%02d' %(loctime.tm_year, loctime.tm_mon, loctime.tm_mday, loctime.tm_hour, loctime.tm_min, loctime.tm_sec)

evtnums = [0,0,0]	# number of events we kept
oevtnums = [0,0,0]	# original number of events
evtdifftimes = [0.,0.,0.]	# alternative to the livetimes entries. Kept just in case
print timestamp
for i in [0]:#range(rqnum):
	firstdone = False
	print i
	rq = {}
	for j in range(len(rqdirstrs[i])):
		# Read one file at a time. This will be slower but should hopefully limit the memory usage
		datapath = os.path.join(rqdir,rqdirstrs[i][j]+'/rootfiles')
		filelist = glob.glob(os.path.join(datapath, '*.rq.root'))
		filelist.sort()
		filenum = len(filelist)
		if max_files_in_dirs[i] != None and filenum > max_files_in_dirs[i]: filenum = max_files_in_dirs[i]
		for m in range(filenum):
			sfile = filelist[m]
			try:
				# unfortunatelly I don't have enough memory to multithread
				froot = rr.read_root(sfile, numberOfThreads=1, max_files_in_dir=max_files_in_dirs[i])
				# read branches
				branches = branches_1d + branches_2d + branches_3d + branches_livetime
				darr = froot.get_branches(branches)
			except:
				print "An exception occured with file %s" %(sfile)
				continue
			# store the original number of events
			if len(darr['pulse_area_phe'].shape) < 2: # single event hence only 1 dimmension. Add one more dimmension
				for branch in branches_1d + branches_2d + branches_3d:
					darr[branch] = n.array([darr[branch]])
			oevtnum = darr['pulse_area_phe'].shape[0]
			if oevtnum == 0: continue	# nothing here
			oevtnums[i] += oevtnum
			# store the difference in event times, just in case.
			evtdifftimes[i] += (darr['luxstamp_samples'].max() - darr['luxstamp_samples'].min())/1e8
			#-------------------- Cuts
			# cut on the luxstamp to remove residual Tritium and DD events
			keep_live_flags = ~ (
				((darr['luxstamp_samples'] > livetime_exclude1[0]) & (darr['luxstamp_samples'] < livetime_exclude1[1])) & 
				((darr['luxstamp_samples'] > livetime_exclude2[0]) & (darr['luxstamp_samples'] < livetime_exclude2[1])) & 
				((darr['luxstamp_samples'] > livetime_exclude3[0]) & (darr['luxstamp_samples'] < livetime_exclude3[1])) & 
				((darr['luxstamp_samples'] > livetime_exclude4[0]) & (darr['luxstamp_samples'] < livetime_exclude4[1]))
				)
			###### Reduce the number of events to store, keep only small enough and
			# large enough events and small events that have only S1s
			# Do our best to store isolated S1s even though the total pulse area is too small
			small_event_s1_flags = ((darr['pulse_classification'] == 1).sum(axis=1) > 0) & (darr['full_evt_area_phe'] < 200)
			# Now keep only the events in the proper area range.
			event_size_flags = (darr['full_evt_area_phe'] > 100) & (darr['full_evt_area_phe'] < 20000)
			# Make sure that these larger events have an S2. We already made sure we
			# captured S1s before. Also do an initial cut of etrains by requiring that
			# there are no more than 4 S2s
			temps2flags = (darr['pulse_classification'] == 2).sum(axis=1)
			event_s2_flags = (temps2flags > 0) & (temps2flags < 5)
			keepflags = keep_live_flags & (small_event_s1_flags | (event_size_flags & event_s2_flags))
			keepinds = n.where(keepflags)[0]
			# store the number of kept events
			evtnum = keepinds.size
			if evtnum == 0: continue	# nothing here
			evtnums[i] += evtnum
			#print evtnum,
			#sys.stdout.flush()
			# store the RQs
			if not firstdone:	# first file, need to initialize the RQ structure
				for branch in branches_livetime:
					rq[branch] = darr[branch]
				for branch in branches_1d+branches_2d:
					rq[branch] = n.take(darr[branch], keepinds, axis=0)
				for branch in branches_3d:
					if branch == 'spike_count':
						tempspikes = n.float32(n.take(darr['spike_count'], botinds, axis=2).sum(axis=2))
						rq['spike_count_bot'] = n.take(tempspikes, keepinds, axis=0)
						tempspikes = n.float32(n.take(darr['spike_count'], topinds, axis=2).sum(axis=2))
						rq['spike_count_top'] = n.take(tempspikes, keepinds, axis=0)
					elif branch == 'peak_area_phe':
						tempspikes = n.take(darr['peak_area_phe'], botinds, axis=2).sum(axis=2)
						rq['pulse_area_phe_bot'] = n.take(tempspikes, keepinds, axis=0)
						tempspikes = n.take(darr['peak_area_phe'], topinds, axis=2).sum(axis=2)
						rq['pulse_area_phe_top'] = n.take(tempspikes, keepinds, axis=0)
					else:
						rq[branch] = n.take(darr[branch].sum(axis=2), keepinds, axis=0)
				firstdone = True
			else: # firstdone
				for branch in branches_livetime:
					rq[branch] = n.append(rq[branch],darr[branch])
				for branch in branches_1d+branches_2d:
					rq[branch] = n.append(rq[branch],n.take(darr[branch], keepinds, axis=0), axis=0)
				for branch in branches_3d:
					if branch == 'spike_count':
						tempspikes = n.float32(n.take(darr['spike_count'], botinds, axis=2).sum(axis=2))
						rq['spike_count_bot'] = n.append(rq['spike_count_bot'],n.take(tempspikes, keepinds, axis=0), axis=0)
						tempspikes = n.float32(n.take(darr['spike_count'], topinds, axis=2).sum(axis=2))
						rq['spike_count_top'] = n.append(rq['spike_count_top'],n.take(tempspikes, keepinds, axis=0), axis=0)
					elif branch == 'peak_area_phe':
						tempspikes = n.take(darr['peak_area_phe'], botinds, axis=2).sum(axis=2)
						rq['pulse_area_phe_bot'] = n.append(rq['pulse_area_phe_bot'],n.take(tempspikes, keepinds, axis=0), axis=0)
						tempspikes = n.take(darr['peak_area_phe'], topinds, axis=2).sum(axis=2)
						rq['pulse_area_phe_top'] = n.append(rq['pulse_area_phe_top'],n.take(tempspikes, keepinds, axis=0), axis=0)
					else:
						rq[branch] = n.append(rq[branch],n.take(darr[branch].sum(axis=2), keepinds, axis=0), axis=0)
			del darr
		print j,
		sys.stdout.flush()
	print ''
	#-------------------------------- Store the data
	# save the arrays making up the RQ. This is a more efficient way to save (less HD usage) and you can select which arrays
	# to read back in
	for key in rq.keys():
		sfile = '%s_backgroundRQs_%s_%s' %(dsets[i], timestamp, key)
		n.save(sfile, rq[key])
	# save auxiulary
	sfile = '%s_backgroundRQs_%s_evtnum-oevtnum-evtdifftime' %(dsets[i], timestamp)
	n.save(sfile, n.array([evtnums[i], oevtnums[i], evtdifftimes[i]]))







