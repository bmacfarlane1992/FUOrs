#
# r1_r.py
#
# Python program to read in rdisc.1 data structures from SEREN analysedisc
# If rdisc files aren't found (for runs with varied Keplerian velocity disc restriction and 
# inclinations), script returns null values.
#
# Author: Benjamin MacFarlane
# Date: 16/03/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
time_check = 0
#
	# Define time refine limits for disc radius/mass plots. As with snaparr, must use dictionary
#
d = {}
tstart2 = 'tstart2' ; tend2 = 'tend2' ; tstart3 = 'tstart3' ; tend3 = 'tend3'
tstart4 = 'tstart4' ; tend4 = 'tend4' ; tstart5 = 'tstart5' ; tend5 = 'tend5'
tstart6 = 'tstart6' ; tend6 = 'tend6' 
#
d[tstart2] = 93.68 ; d[tend2] = 94.58 ; d[tstart3] = 88.00 ; d[tend3] = 94.00
d[tstart4] = 93.93 ; d[tend4] = 94.93 ; d[tstart5] = 91.50 ; d[tend5] = 93.50
d[tstart6] = 86.00 ; d[tend6] = 89.00 
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
import os
import random
import math
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.axes as axes
import pylab
from astropy.io import ascii
import glob
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
def read(arch_dir, plotdir, ea_run, snaparr, v_K, inclin):
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   # Read rdisc.1 file if original DS runs are being analysed #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	if ((v_K == "90") and (inclin == "0")):
		print "rdisc.1 file being read/plots being generated"
		arch_ext = 'rdisc_DS/DE05.rdisc.1'
		filename = arch_dir+arch_ext
#
	# Define arrays to fill
#
		time = [] ; m_s = [] ; m_mri_d = [] ; r_d_kep = [] ; m_d_kep = [] ; r_d_piv = []
		m_d_piv = [] ; r_d_sigALMA = [] ; m_d_sigALMA = [] ; snap = [] ; dtheta = []
		mcloud1 = [] ; mcloud5 = [] ; mcloud10 = []
#
	# Extra parameters for ea runs with episodic accretion
#
		acctag_num = [] ; acctag = []
#
	# Read in ea runs with episodic accretion
#
		if (float(ea_run) > 1):
#
			f = open(filename, 'r')
			for line in f:
				line = line.strip()
				columns = line.split()
				time.append(float(columns[0])/1000.) ; m_s.append(float(columns[1]))
				m_mri_d.append(float(columns[2])) ; r_d_kep.append(float(columns[3]))
				m_d_kep.append(float(columns[4])) ; r_d_piv.append(float(columns[5]))
				m_d_piv.append(float(columns[6])) ; r_d_sigALMA.append(float(columns[11]))
				m_d_sigALMA.append(float(columns[12])) ; snap.append(columns[13])
				acctag_num.append(float(columns[14])) ; acctag.append(columns[15])
				dtheta.append(float(columns[16])) ; mcloud1.append(float(columns[17]))
				mcloud5.append(float(columns[18])) ; mcloud10.append(float(columns[19]))
			f.close()
#
	# Read in ea runs without episodic accretion
#
		if (float(ea_run) <= 1):
#
			f = open(filename, 'r')
			for line in f:
				line = line.strip()
				columns = line.split()
				time.append(float(columns[0])/1000.) ; m_s.append(float(columns[1]))
				m_mri_d.append(float(columns[2])) ; r_d_kep.append(float(columns[3]))
				m_d_kep.append(float(columns[4])) ; r_d_piv.append(float(columns[5]))
				m_d_piv.append(float(columns[6])) ; r_d_sigALMA.append(float(columns[11]))
				m_d_sigALMA.append(float(columns[12])) ; snap.append(columns[13])
				dtheta.append(float(columns[14])) ; mcloud1.append(float(columns[15]))
				mcloud5.append(float(columns[16])) ; mcloud10.append(float(columns[17]))
#
				acctag_num.append(0)
#
			f.close()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# Carry out analysis of file in preparation for plotting #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
#
	# Determine the number of times at which accretion is occuring in prep. for overplotting
	# lines onto mass and radius evolution
#
		acc_count = 0
		for i in range(0, int(len(time))):
			if (int(acctag_num[i]) == 1):
				acc_count = acc_count + 1
#
	# Calculate ABSOLUTE differences in mass and size values estimated from either Keplerian
	# or azimuthal criteria, by constructing array of size n_snaps - number of time snapshots
#
		n_snaps = int(len(time))
		m_diff_kphi = [0.] * n_snaps ; m_diff_ksig = [0.] * n_snaps
		m_diff_phisig = [0.] * n_snaps ; r_diff_kphi = [0.] * n_snaps
		r_diff_ksig = [0.] * n_snaps ; r_diff_phisig = [0.] * n_snaps
#
		for i in range(0,n_snaps):
			m_diff_kphi[i] = abs(m_d_kep[i] - m_d_piv[i])
			m_diff_ksig[i] = abs(m_d_kep[i] - m_d_sigALMA[i])
			m_diff_phisig[i] = abs(m_d_piv[i] - m_d_sigALMA[i])
			r_diff_kphi[i] = abs(r_d_kep[i] - r_d_piv[i])
			r_diff_ksig[i] = abs(r_d_kep[i] - r_d_sigALMA[i])
			r_diff_phisig[i] = abs(r_d_piv[i] - r_d_sigALMA[i])
#
	# Print time limits for individual runs for refining of times plotted
#
		if (time_check == 1):
			print "Minimum time of run is: ", min(time), " kyr"
			print "Maximum time of run is: ", max(time), " kyr"
			exit()
#
	# Define time indices for accretion times, in order to generate continuous fill
#
		hasharr = []
		jstart = 0
		for i in range(0, n_snaps-1):
			for j in range(jstart, n_snaps-1):
				if (int(acctag_num[j]) == 1):
					hasharr.append(j)
					for k in range(j, n_snaps):
						if(int(acctag_num[k]) != int(acctag_num[j])):
							hasharr.append(k) ; jstart = k ; break
					break	
		hasharr_app=[]
		for i in hasharr:
 			if i not in hasharr_app:
       				hasharr_app.append(i)
		if (acctag_num[n_snaps-1] == "1"):
			hasharr_app.append(n_snaps-1)	
		n_accr = len(hasharr_app)/2
		hasharr_app = np.reshape(hasharr_app, (n_accr, 2))
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# Plot rdisc.1 temporally evolved parameters if file is complete #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# Plot mass of star and MRI disc vs. time
#
		plt.figure(1)
		plt.subplot(211)
		line1 = plt.plot(time, m_s)
		plt.ylabel("Stellar mass "+(r'(M$_{\odot}$)'))
		plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
#	
		plt.subplot(212)
		line1 = plt.plot(time, m_mri_d)
		plt.xlabel('time (kyr)')
		pylab.xlim([int(min(time)),int(max(time))])
		plt.ylabel("MRI disc mass "+(r'(M$_{\odot}$)'))
#	
		plt.savefig(plotdir+'iad_star_mass.pdf')
		plt.clf()
#
	# Plot mass of disc with Keplerian and azimuthal velocity criteria
#
		plt.figure(1)
		ax1 = plt.subplot(211)
		line11 = plt.plot(time, m_d_kep, label="Keplerian")
		line22 = plt.plot(time, m_d_piv, label="Azimuthal")
		line33 = plt.plot(time, m_d_sigALMA, label="ALMA SD")
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			0, ax1.get_ylim()[1], color='k', alpha = 0.5)
		legend = plt.legend(loc = 'upper left', fontsize=8)
		plt.xlabel('time (kyr)')
		pylab.xlim([int(min(time)),math.ceil(max(time))])
		plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1])
		plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
		plt.ylabel("disc mass "+(r'(M$_{\odot}$)'))
#
	# Plot absolute difference of disc mass when different criteria are applied
#
		ax2 = plt.subplot(212)
		line1 = plt.plot(time, m_diff_kphi, label="Keplerian vs. Azimuthal")
		line2 = plt.plot(time, m_diff_ksig, label="Keplerian vs. ALMA SD")
		line3 = plt.plot(time, m_diff_phisig, label="Azimuthal vs. ALMA SD")
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			0, ax2.get_ylim()[1], color='k', alpha = 0.5)
		legend = plt.legend(loc = 'upper left', fontsize=8)
		plt.xlabel('time (kyr)')
		pylab.xlim([int(min(time)),math.ceil(max(time))])
		plt.ylim(0, ax2.get_ylim()[1])
		plt.ylabel(('|$\Delta$ ')+"M| (AU)")
		plt.savefig(plotdir+'disc_mass.pdf')	
		plt.clf()
#
	# Plot radius of disc with Keplerian and azimuthal velcity criteria
#
		plt.figure(1)
		ax1 = plt.subplot(211)
		line11 = plt.plot(time, r_d_kep, label="Keplerian")
		line22 = plt.plot(time, r_d_piv, label="Azimuthal")
		line33 = plt.plot(time, r_d_sigALMA, label="ALMA SD")
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			0, ax1.get_ylim()[1], color='k', alpha = 0.5)
		legend = plt.legend(loc = 'upper left', fontsize=8)
		plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
		pylab.xlim([int(min(time)), math.ceil(max(time))])
		plt.ylim(0, ax1.get_ylim()[1])
		plt.ylabel("disc radius (AU)")
#
	# Plot absolute difference of disc radius when different criteria are applied
#
		ax2 = plt.subplot(212)
		line1 = plt.plot(time, r_diff_kphi, label="Keplerian vs. Azimuthal")
		line2 = plt.plot(time, r_diff_ksig, label="Keplerian vs. ALMA SD")
		line3 = plt.plot(time, r_diff_phisig, label="Azimuthal vs. ALMA SD")
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			0, ax2.get_ylim()[1], color='k', alpha = 0.5)
		legend = plt.legend(loc = 'upper left', fontsize=8)
		plt.xlabel('time (kyr)')
		pylab.xlim([int(min(time)),math.ceil(max(time))])
		plt.ylim(0, max(r_diff_kphi))
		plt.ylabel(('|$\Delta$ ')+"R| (AU)")
		plt.savefig(plotdir+'disc_radius.pdf')	
		plt.clf()
#
	# Plot refined time slice, only for EA runs to evaluate disc parameters
	# around an outburst event
#
		if (snaparr.ndim == 2):
#
			tstart = d['tstart'+str(ea_run)]
			tend = d['tend'+str(ea_run)]
#
			snapstart = int((tstart - min(time)) / ((max(time) - min(time)) / n_snaps))
			snapend = int((tend - min(time)) / ((max(time) - min(time)) / n_snaps))
#
	# Plot time restricted mass of disc with Keplerian and azimuthal velocity criteria
#
			plt.figure(1)
			ax1 = plt.subplot(111)
			line11 = plt.plot(time[snapstart:snapend], m_d_kep[snapstart:snapend], \
			   label="Keplerian")
			line22 = plt.plot(time[snapstart:snapend], m_d_piv[snapstart:snapend], \
			   label="Azimuthal")
			line33 = plt.plot(time[snapstart:snapend], m_d_sigALMA[snapstart:snapend], \
			   label="ALMA SD")
			for i in range(0,n_accr):
				plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
				0, ax1.get_ylim()[1], color='k', alpha = 0.5)
			legend = plt.legend(loc = 'upper left', fontsize=8)
			plt.xlabel('time (kyr)')
			pylab.xlim([min(time[snapstart:snapend]),max(time[snapstart:snapend])])
			plt.ylim(0, ax1.get_ylim()[1])
			plt.ylabel("disc mass "+(r'(M$_{\odot}$)'))
			plt.savefig(plotdir+'trestr_disc_mass.pdf')	
			plt.clf()
#
	# Plot time restricted radius of disc with Keplerian and azimuthal velocity criteria
#
			plt.figure(1)
			ax1 = plt.subplot(111)
			line11 = plt.plot(time[snapstart:snapend], r_d_kep[snapstart:snapend], \
			   label="Keplerian")
			line22 = plt.plot(time[snapstart:snapend], r_d_piv[snapstart:snapend], \
			   label="Azimuthal")
			line33 = plt.plot(time[snapstart:snapend], r_d_sigALMA[snapstart:snapend], \
			   label="ALMA SD")
			for i in range(0,n_accr):
				plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
				0, ax1.get_ylim()[1], color='k', alpha = 0.5)
			legend = plt.legend(loc = 'upper left', fontsize=8)
			plt.xlabel('time (kyr)')
			pylab.xlim([min(time[snapstart:snapend]),max(time[snapstart:snapend])])
			plt.ylim(0, ax1.get_ylim()[1])
			plt.ylabel("disc radius (AU)")
			plt.savefig(plotdir+'trestr_disc_radius.pdf')	
			plt.clf()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   # Write hasharr_app data to ensure mass_comp script successfully traces accretion events #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
		f = open(arch_dir+'../acc_params.dat','w')
		for i in range(0, n_accr):
			f.write( str(hasharr_app[i][0]*0.01 + min(time))+' '+str(hasharr_app[i][1]*0.01 + min(time))+'\n' )
		f.close()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   # Return null values for non-disc mass criterion values, if rdisc.1 file incomplete#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	else:
		print "rdisc.1 file not present in this run, returning null values"
		acc_count = 0 ; acctag_num = 0 ; time = 0 ; hasharr_app = 0 ; n_accr = 0
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   # Read in rdisc parameters of disc mass based on v_phi, orbital dynamic and density criteria #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	print "Mass criteria for snapshots now being analysed"
	arch_ext = 'rdisc/DE05.rdisc.1'
	filename = arch_dir+arch_ext
#
	# Define arrays to fill, read file and store data in relevant array
#
	m_s = [] ; m_mri_d = [] ; m_d_kep = [] ; m_d_piv = [] ; m_d_sigALMA = [] ; r_d_kep = []
#
	f = open(filename, 'r')
	for line in f:
		line = line.strip()
		columns = line.split()
		m_s.append(float(columns[1])) ; m_mri_d.append(float(columns[2])) 
		r_d_kep.append(float(columns[3])) ; m_d_kep.append(float(columns[4]))
		m_d_piv.append(float(columns[6])) ; m_d_sigALMA.append(float(columns[12]))
	f.close()
	
#	
	return acc_count, acctag_num, hasharr_app, n_accr, r_d_kep, \
	   m_s, m_mri_d, m_d_kep, m_d_piv, m_d_sigALMA
#
#
