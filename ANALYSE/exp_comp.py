#
# exp_comp.py
#
# Programme to plot comparison of radial T and SD exponents for no, continuous and episodic
# feedback regimes
#
# Author: Benjamin MacFarlane
# Date: 18/04/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
snaparr = [[131,331],[600, 1100, 2100],[1050, 1150, 1350, 1450]]		# Select snapshots (from respective snaparr indices) 
ea_run = ["NF", "CF", "EF"]							# from which to compare exponents from
#										# Also select which EA run to draw exponents
#  										# from (see main.py for selection available)
#
col_arr = ['b','g','r','c']
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
def comp(arch_dir, plotdir):
#
	print "\n Exponent comparisons ([Sigma, T] profiles) between feedback regimes is now being plotted \n"
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Read in files with exponent data
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	SD_ea = [] ; SD_snap = [] ; SD_exp = []
	f = open(arch_dir+'SD_exps.dat','r')
	for line in f:
		line = line.strip() ; columns = line.split()
		SD_ea.append(float(columns[0])) ; SD_snap.append(float(columns[1]))
		SD_exp.append(float(columns[2]))
	f.close()
	SD_ea = np.array(SD_ea) ; SD_snap = np.array(SD_snap) ; SD_exp = np.array(SD_exp)
#
	T_ea = [] ; T_snap = [] ; T_exp = []
	f = open(arch_dir+'T_exps.dat','r')
	for line in f:
		line = line.strip() ; columns = line.split()
		T_ea.append(float(columns[0])) ; T_snap.append(float(columns[1]))
		T_exp.append(float(columns[2]))
	f.close()
	T_ea = np.array(T_ea) ; T_snap = np.array(T_snap) ; T_exp = np.array(T_exp)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Plotting of exponent comparisons for different feedback regimes
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	# SD exponent comparison
#
	plt.figure(1)
#
	ax1 = plt.subplot(111)
#
	for i in range(0,len(ea_run)):
		for j in range(0,len(snaparr[i])):
			for k in range(0,len(SD_ea)):
				if (snaparr[i][j] == SD_snap[k]):
					plt.scatter((i),abs(SD_exp[k]), s=80, facecolors=col_arr[j], edgecolors=col_arr[j])
					plt.xticks(range(len(ea_run)), ea_run)
	ax1.add_patch( \
	   patches.Rectangle((ax1.get_xlim()[0], 1.), \
           (ax1.get_xlim()[1] - ax1.get_xlim()[0]), 0.5, \
	   color='k', alpha = 0.5 ) )
	plt.xlabel('Feedback regime' ) ; plt.ylabel('p')
#
	plt.savefig(plotdir+'SD_exp.pdf') ; plt.clf()
#
	# T exponent comparison
#
	plt.figure(1)
#
	ax1 = plt.subplot(111)
#
	for i in range(0,len(ea_run)):
		for j in range(0,len(snaparr[i])):
			for k in range(0,len(T_ea)):
				if (snaparr[i][j] == T_snap[k]):
					plt.scatter((i),abs(T_exp[k]), s=80, facecolors=col_arr[j], edgecolors=col_arr[j])
					plt.xticks(range(len(ea_run)), ea_run)
	ax1.add_patch( \
	   patches.Rectangle((ax1.get_xlim()[0], 0.35), \
           (ax1.get_xlim()[1] - ax1.get_xlim()[0]), 0.45, \
	   color='k', alpha = 0.5 ) )
	plt.xlabel('Feedback regime' ) ; plt.ylabel('q')
#
	plt.savefig(plotdir+'T_exp.pdf') ; plt.clf()

