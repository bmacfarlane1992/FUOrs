#
# rdisc_r.py
#
# Python program to read in rdisc files, and plot temporally evolved parameters
#
# Author: Benjamin MacFarlane
# Date: 28/01/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
import random
import math
import numpy as np
import sys
import matplotlib.pyplot as plt
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
def read(arch_dir, plotdir, acc_run, hasharr_app, n_accr, r_inspec, v_K, inclin): 	
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Read rdisc.1.i files, if original DS run being analysed #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	if ((v_K == "90") and (inclin == "0")):
		print("rdisc files being read")
		filename = arch_dir+'rdisc_DS/DE05.rdisc.1.'+str(r_inspec)
#
	# Define arrays to fill
#
		time = [] ; r = [] ; Q = [] ; T = [] ; sig = [] ; Omeg = [] ; DiscM = []
		rInM = [] ; M_s = [] ; M_d = [] ; v_r = [] ; v_the = [] ; v_z = [] ; vkep = []
		h = [] ; h_mid = [] ; rpart = []
#
	# Loop over files in pdisc directory to read all temporally evolved values
	# Ignore header of file whilst reading
#
		f = open(filename, 'r')
		header = f.readline()
		for line in f:
			line = line.strip()
			columns = line.split()
			time.append(float(columns[0])/1000.) ;	r.append(float(columns[1]))
			Q.append(float(columns[2])) ; T.append(float(columns[3]))
			sig.append(float(columns[4])) ; Omeg.append(float(columns[5]))
			DiscM.append(float(columns[6])) ; rInM.append(float(columns[7]))
			M_s.append(float(columns[8])) ; M_d.append(float(columns[9]))
			v_r.append(float(columns[10])) ; v_the.append(float(columns[11]))
			v_z.append(float(columns[12])) ; vkep.append(float(columns[13]))
			h.append(float(columns[14])) ; h_mid.append(float(columns[15]))
			rpart.append(float(columns[16]))
		f.close()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Plotting temporally evolved disc parameters for set radius
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
		print("Plotting rdisc parameters for "+str(r_inspec)+" AU")
#
		fig = plt.figure(1, figsize=(60,40))
		fig.subplots_adjust(hspace=.3, wspace = 0.4)
		fig.suptitle('Disc parameters at '+str(r_inspec)+' AU', fontsize=12)
#
		ax1 = plt.subplot(231)
		line1 = plt.plot(time, Q)
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			0, ax1.get_ylim()[1], color='k', alpha = 0.5)
		plt.xlabel("Time (kyr)")
		pylab.xlim([int(min(time)),int(max(time))])
		if (max(Q) < 20):
			plt.ylim(0, max(Q))
		else:
			plt.ylim(0, 20)
		plt.ylabel("Toomre Q parameter")
		plt.yticks(fontsize = 8)
		plt.xticks(fontsize = 8)
#
		ax2 = plt.subplot(232)
		line2 = plt.plot(time, T)
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			0, max(T), color='k', alpha = 0.5)
		plt.xlabel("Time (kyr)")
		pylab.xlim([int(min(time)),int(max(time))])
		plt.ylim(0, max(T))	
		plt.ylabel("Temperature (K)")
		plt.yticks(fontsize = 8)
		plt.xticks(fontsize = 8)
#
		ax3 = plt.subplot(233)
		line3 = plt.plot(time, sig)
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			0, max(sig), color='k', alpha = 0.5)
		plt.xlabel("Time (kyr)")
		pylab.xlim([int(min(time)),int(max(time))])
		plt.ylim(0, max(sig))	
		plt.ylabel("Surface density")
		plt.yticks(fontsize = 8)
		plt.xticks(fontsize = 8)
#
		ax4 = plt.subplot(234)
		line4 = plt.plot(time, DiscM)
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			0, ax4.get_ylim()[1], color='k', alpha = 0.5)
		plt.xlabel("Time (kyr)")
		pylab.xlim([int(min(time)),int(max(time))])
		plt.ylim(0, ax4.get_ylim()[1])
		plt.ylabel("Disc mass "+(r'(M$_{\odot}$)'))
		plt.yticks(fontsize = 8)
		plt.xticks(fontsize = 8)
#
		ax5 = plt.subplot(235)
		line5 = plt.plot(time, v_the)
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			ax5.get_ylim()[0], ax5.get_ylim()[1], color='k', alpha = 0.5)
		plt.xlabel("Time (kyr)")	
		pylab.xlim([int(min(time)),int(max(time))])
		plt.ylim(ax5.get_ylim()[0], ax5.get_ylim()[1])	
		plt.ylabel((r'v$_{\phi}$')+' (km'+(r's$^{-1}$')+')')
		plt.yticks(fontsize = 8)
		plt.xticks(fontsize = 8)
#
		ax6 = plt.subplot(236)
		line6 = plt.plot(time, v_r)
		for i in range(0,n_accr):
			plt.fill_between(time[hasharr_app[i][0]:hasharr_app[i][1]], \
			ax6.get_ylim()[0], ax6.get_ylim()[1], color='k', alpha = 0.5)
		plt.xlabel("Time (kyr)")
		pylab.xlim([int(min(time)),int(max(time))])
		plt.ylim(ax6.get_ylim()[0], ax6.get_ylim()[1])	
		plt.ylabel((r'v$_{r}$')+' (km'+(r's$^{-1}$')+')')
		plt.yticks(fontsize = 8)
		plt.xticks(fontsize = 8)
#
		plt.savefig(plotdir+str(r_inspec)+'AU_evol.png')
		plt.clf()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# If EA run is not from original data, return to main program#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	else:
		print "rdisc files not present in this run, returning to main.py"
#
#
