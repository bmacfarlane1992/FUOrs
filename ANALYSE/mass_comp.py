#
# mass_comp.py
#
# Programme to plot comparison of masses in simulation and PV diagram analyses
#
# Author: Benjamin MacFarlane
# Date: 19/02/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
print_term = "FALSE"
#
t_ref = 2.5
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
import matplotlib.pyplot as cm
import scipy
from scipy.optimize import curve_fit
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
def comp(arch_dir, plotdir, ea_run, snaparr, timearr, v_K, inclin, m_s, pmass, m_d_kep, \
   m_d_piv, m_d_sigALMA, pv_mass):
#
	print "Mass comparisons of simulation and PV data now being plotted"
#
	# Define number of files to be read dependent on snaparr dimensions
	# Set snaparr reference to 1D using snaparr_tmp dummy
	# Also, restructure pmass array into pmass_tmp dummy, for consistency to file_n pointer.
	# pmass_tmp stores sum of planetary mass to evaluate star + planet + disc mass properly
#
	# For EA [0, 1] runs
#
	if (snaparr.ndim == 1):
		file_n = len(snaparr)
		snaparr_tmp = [0]*file_n
		timearr_tmp = [0]*file_n
		pmass_tmp = [0]*file_n
		fcount = 0
		for i in range(0, len(snaparr)):
			snaparr_tmp[fcount] = snaparr[i]
			timearr_tmp[fcount] = timearr[i]
			for a in range(0, len(pmass)):
				pmass_tmp[fcount] = pmass_tmp[fcount] + pmass[a][i]
			fcount = fcount + 1
#
	# For EA [2, 3, 4, 5, 6] runs
#
	elif (snaparr.ndim == 2):
		file_n = len(snaparr)*len(snaparr[0])
		snaparr_tmp = [0]*file_n
		timearr_tmp = [0]*file_n
		pmass_tmp = [0]*file_n
		fcount = 0
		for i in range(0, len(snaparr)):
			for j in range(0, len(snaparr[0])):
				snaparr_tmp[fcount] = snaparr[i][j]
				timearr_tmp[fcount] = timearr[i][j]
				for a in range(0, len(pmass)):
					pmass_tmp[fcount] = pmass_tmp[fcount] + pmass[a][i][j]
				fcount = fcount + 1
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Loop over snapshots, and compute masses for comparison #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	m_sys_kep = [0] * file_n
	m_sys_piv = [0] * file_n
	m_sys_sigALMA = [0] * file_n
#
	for i in range(0, file_n):
#
		m_sys_kep[i] = m_s[i] + pmass_tmp[i] + m_d_kep[i]
		m_sys_piv[i] = m_s[i] + pmass_tmp[i] + m_d_piv[i]
		m_sys_sigALMA[i] = m_s[i] + pmass_tmp[i] + m_d_sigALMA[i]
#
		if (print_term == "TRUE"):
			print "for snaparr value: "+str(snaparr_tmp[i])
			print "\tMass of central protostar: "+str(m_s[i])
			print "\tMass of planets in disc: "+str(pmass_tmp[i])
			print "\tMass of system (Keplerian velocity criterion): "+str(m_sys_kep[i])
			print "\tMass of system (Oribit infall criterion): "+str(m_sys_piv[i])
			print "\tMass of system (ALMA density criterion): "+str(m_sys_sigALMA[i])
			print "\tMass of system (Keplerian fitted P-V diagram): "+str(pv_mass[i])
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Define simulation independent time dimension, tau to evaluate differences between #
	# simulation and PV determined system masses
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	t_len = [] ; t_cen = [] ; t_b1 = [] ; t_b2 = [] ; t_s = [] ; t_e = []
#
	if (snaparr.ndim == 2):
#
		f = open('../DATA/'+str(ea_run)+'/acc_params.dat', 'r')
		for line in f:
			line = line.strip() ; columns = line.split()
			t_len.append(float(columns[0])) ; t_cen.append(float(columns[1]))
			t_b1.append(float(columns[2])) ; t_b2.append(float(columns[3]))
			t_s.append(float(columns[4])) ; t_e.append(float(columns[5]))
		f.close()

#
		tau = [0] * file_n
		for i in range(0, file_n):
			for j in range(0, len(t_len)):	
				if ( (timearr_tmp[i] > t_b1[j]) and \
				   (timearr_tmp[i] < t_b2[j]) ):
#
					tau[i] = ( (timearr_tmp[i] + t_len[j] - t_ref) / t_cen[j]) - 1.
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Plotting of mass comparisons, for different accretion events
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
		plt.figure(1)
		ax1 = plt.subplot(111)
#
		plt.plot(timearr_tmp, pv_mass, linestyle = 'dashed', color = 'b', label = 'PV')
		plt.plot(timearr_tmp, m_sys_kep, linestyle = 'solid', color = 'g', label = 'Keplerian')
		plt.plot(timearr_tmp, m_sys_piv, linestyle = 'solid', color = 'r', label = 'Orbit infall')
		plt.plot(timearr_tmp, m_sys_sigALMA, linestyle = 'solid', color = 'k', label = 'ALMA SD')
		plt.ylim(0, ax1.get_ylim()[1])
		for i in range(0, len(t_len)):
			plt.fill_between((t_s[i], t_e[i]), \
			   0, ax1.get_ylim()[1], color='k', alpha = 0.5)
#
		legend = plt.legend(loc = 'upper left', fontsize=8)
		plt.xlabel('Time, '+(r'$\tau$') )
		plt.ylabel('Mass '+(r'(M$_{\odot}$)') )
		plt.savefig(plotdir+'mass_comp.png')
		plt.clf()
#
#
