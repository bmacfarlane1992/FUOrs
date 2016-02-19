#
# sink_r.py
#
# Python program to read in DE05.sink[i] data structures from SEREN analysedisc
# executable. Uses sink data from disc to output sink parameters for plotting routines in
# both pdisc_EAr and pdisc_r routines
#
# Author: Benjamin MacFarlane
# Date: 29/01/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
pcAU = 206265.
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
def read(arch_dir, plotdir, ea_run, snaparr):
	print "Planet (sink) data being read"
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Read sink files, defining arrays to be filled for sink parameters #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	pmass = []
	pradius = []
	ptime = []
	file_list = sorted(glob.glob(arch_dir+'../../../ICs/'+str(ea_run)+'/DE05.sink*'))
	file_n = len(file_list)
#
	# Define arrays to fill
#
	time = [[] for i in range(file_n)]
	mass = [[] for i in range(file_n)]
	radius = [[] for i in range(file_n)]
#
	for i in range(0,file_n):		
#
	# Now loop over lines in each file and extract variables of interest
	# then convert to numpy format for future data extraction for plotting
#
		if (snaparr.ndim == 2):
#
			f = open(file_list[i], 'r')
			for line in f:
				line = line.strip() ; columns = line.split()
				time[i].append(float(columns[1])*1000.)
				mass[i].append(float(columns[8]))
				radius[i].append(math.sqrt((float(columns[2])*pcAU/10.)**2 + \
			   	   (float(columns[3])*pcAU/10.)**2 + \
				   (float(columns[4])*pcAU/10.)**2))
			f.close()
#
		elif (snaparr.ndim == 1):
#
			f = open(file_list[i], 'r')
			for line in f:
				line = line.strip() ; columns = line.split()
				time[i].append(float(columns[1])*1000.)
				mass[i].append(float(columns[8]))
				radius[i].append(math.sqrt((float(columns[2])*pcAU/10.)**2 + \
   				   (float(columns[3])*pcAU/10.)**2 + \
				   (float(columns[4])*pcAU/10.)**2))
			f.close()

	time = np.array(time) ; mass = np.array(mass) ; radius = np.array(radius)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Extract sink radius/mass for snapshot selection #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	timearr = []
#
	# For EA runs [0, 1]
#
	if (snaparr.ndim == 1):
#
		for i in range(len(snaparr)):
			timearr.append( (snaparr[i] * 0.01) + float(min(time[0])) )
		timearr = np.array(timearr)
#
		pmass = [ [ 0 for j in xrange(len(timearr)) ] for i in xrange(len(time)) ]
		pradius = [ [ 0 for j in xrange(len(timearr)) ] for i in xrange(len(time)) ]
		ptime = [ [ 0 for j in xrange(len(timearr)) ] for i in xrange(len(time)) ]
		for a in range(0,len(time)):	# Loop over number of sinks
			for i in range(0,len(timearr)):	# Loop over time of snapshot under analysis
				for t_s in range(0, len(time[a])):	# Loopover timesnaps in sink file
					if (round(timearr[i], 3) == round(time[a][t_s], 3) ):
						pmass[a][i] = mass[a][t_s]
						pradius[a][i] = radius[a][t_s]
						ptime[a][i]  = time[a][t_s]
#
	# For EA runs [2, 3, 4, 5, 6]
#
	if (snaparr.ndim == 2):
#
		for i in range(0,len(snaparr)):
			for j in range(0,len(snaparr[i])):
				timearr.append( (snaparr[i][j] * .01) + float(min(time[0]))  )
		timearr = np.array(timearr)
		timearr = np.reshape(timearr, (len(snaparr),len(snaparr[i])) )
#
		pmass = [[[0 for k in xrange(len(timearr[0]))] \
		   for j in xrange(len(timearr))] for i in xrange(len(time))]
		pradius = [[[0 for k in xrange(len(timearr[0]))] \
		   for j in xrange(len(timearr))] for i in xrange(len(time))]
		ptime = [[[0 for k in xrange(len(timearr[0]))] \
		   for j in xrange(len(timearr))] for i in xrange(len(time))]
		for a in range(0,len(time)):							    # Loop over number of sinks
			for i in range(0,len(timearr[0])):					    # Loop over accretion events under analysis
				for j in range(0,len(timearr)): 				    # Loop over reference time to accr. event
					for t_s in range(0, len(time[a])):			    # Loop over timesnaps in sink file
						if (round(timearr[i][j],3)==round(time[a][t_s],3)):
							pmass[a][i][j] = mass[a][t_s]
							pradius[a][i][j] = radius[a][t_s]
							ptime[a][i][j]  = time[a][t_s]
#
#
	return pmass, pradius, timearr
#
#
