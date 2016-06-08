#
# exp_comp.py
#
# Programme to plot comparison of radial T and SD exponents for no, continuous and episodic
# feedback regimes
#
# Author: Benjamin MacFarlane
# Date: 08/06/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# Select snapshots (from respective snaparr indices) from which to compare exponents from.  
	# Also select which EA run to draw exponents from (see main.py for selection available)

snaparr = [[310, 360, 410, 460, 510, 560, 610, 660, 710, 760, 780, 800], \
	[150, 300, 450, 600, 750, 900, 1050, 1200, 1350, 1500, 1650, 1800], \
	[100, 150, 200, 250, 400, 450, 550, 600, 1050, 1150, 1350, 1450], \
	[200, 209, 212, 220, 860, 875, 880, 900, 1525, 1585, 1595, 1625], \
	[150, 183, 191, 224, 719, 755, 765, 801, 1342, 1377, 1399, 1434] ]
ea_run = ["NF", "CF", "EF-A","EF-B","EF-C"]		
ea_run_n = [0, 1, 3, 4, 5]
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
def comp(dat_dir, plotdir, col_arr):
#
	print "\n Exponent comparisons ([Sigma, T] profiles) between feedback regimes is now being plotted \n"
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Read in files with exponent data
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	SD_ea = [] ; SD_snap = [] ; SD_time = [] ; SD_exp = [] ; exp_err = []
	f = open(dat_dir+'SD_exps.dat','r')
	for line in f:
		line = line.strip() ; columns = line.split()
		SD_ea.append(float(columns[0])) ; SD_snap.append(float(columns[1]))
		SD_time.append(float(columns[2])) ; SD_exp.append(float(columns[3]))
		exp_err.append(float(columns[4]))
	f.close()
	SD_ea = np.array(SD_ea) ; SD_snap = np.array(SD_snap) ; SD_time = np.array(SD_time)
	SD_exp = np.array(SD_exp) ; exp_err = np.array(exp_err)
#
	T_ea = [] ; T_snap = [] ; T_time = [] ; T_exp = []
	f = open(dat_dir+'T_exps.dat','r')
	for line in f:
		line = line.strip() ; columns = line.split()
		T_ea.append(float(columns[0])) ; T_snap.append(float(columns[1]))
		T_time.append(float(columns[2])) ; T_exp.append(float(columns[3]))
	f.close()
	T_ea = np.array(T_ea) ; T_snap = np.array(T_snap)
	T_time = np.array(T_time) ; T_exp = np.array(T_exp)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Plotting of exponent comparisons for different feedback regimes
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	# SD exponent 
#
	fig = plt.figure(1)
	axcount = 0
	ytck = [0., 1., 2., 3., 4.]
#
	for i in range(0,len(ea_run)):
#
		axcount = axcount+1
		axpoint = 510 + axcount
		ax1 = plt.subplot(axpoint)
#
		if (i != len(ea_run)-1):
			ax1.xaxis.set_major_formatter(plt.NullFormatter())
#
	# Read in accretion parameters and specify T_snap and SD_snap indices are in outburst
#
		n_accr = 0
		if (ea_run_n[i] >  1):
			t_s = [] ; t_e = []
			f = open(dat_dir+str(ea_run_n[i])+'/acc_params.dat','r')
			for line in f:
				line = line.strip() ; columns = line.split()
				t_s.append(float(columns[0])) ; t_e.append(float(columns[1]))
			f.close()
			t_s = np.array(t_s) ; t_e = np.array(t_e) ; n_accr = len(t_s)
#
		n_out = 0 ; y_out = 0
#
		for j in range(0,len(snaparr[i])):
			for k in range(0,len(SD_ea)):
#
	# Define which EA snapshots occur during outburst phase
#
				outburst = 0
				for a in range(0, n_accr):
					if ((SD_ea[k] > 1) and (SD_time[k] > t_s[a]) and (SD_time[k] < t_e[a])):
						outburst = 1
#
				if ((snaparr[i][j] == SD_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 0) and (outburst == 0)):
					plt.scatter(SD_time[k], abs(SD_exp[k]), s=80, \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
					ax1.set_yticks(ytck)
					n_out = n_out + 1
				if ((snaparr[i][j] == SD_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 1) and (outburst == 0)):
					plt.scatter(SD_time[k], abs(SD_exp[k]), s=80, marker="s", \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
					ax1.set_yticks(ytck)
					n_out = n_out + 1
				if ((snaparr[i][j] == SD_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 0) and (outburst == 1)):
					plt.scatter(SD_time[k], abs(SD_exp[k]), s=80, \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
					ax1.set_yticks(ytck)
					y_out = y_out + 1
				if ((snaparr[i][j] == SD_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 1) and (outburst == 1)):
					plt.scatter(SD_time[k], abs(SD_exp[k]), s=80, marker="s", \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
					ax1.set_yticks(ytck)
					y_out = y_out + 1
				plt.ylim(0, 4.) ; plt.xlim(78., 100.) ; plt.ylabel('p')
			plt.legend(loc='upper right', fontsize = 10, scatterpoints = 1)
	plt.xlabel('Time (kyr)' )
#
	plt.savefig(plotdir+'SD_exp.pdf') ; plt.clf()
#
	# T exponent comparison
#
	plt.figure(1)
	axcount = 0
	ytck = [0, 0.75, 1.5, 2.25]
#
	for i in range(0,len(ea_run)):
#
		axcount = axcount+1
		axpoint = 510 + axcount
		ax1 = plt.subplot(axpoint)
#
		if (i != len(ea_run)-1 ):
			ax1.xaxis.set_major_formatter(plt.NullFormatter())
#
	# Read in accretion parameters and specify T_snap and SD_snap indices are in outburst
#
		n_accr = 0
		if (ea_run_n[i] >  1):
			t_s = [] ; t_e = []
			f = open(dat_dir+str(ea_run_n[i])+'/acc_params.dat','r')
			for line in f:
				line = line.strip() ; columns = line.split()
				t_s.append(float(columns[0])) ; t_e.append(float(columns[1]))
			f.close()
			t_s = np.array(t_s) ; t_e = np.array(t_e) ; n_accr = len(t_s)
#
		n_out = 0 ; y_out = 0
#
		for j in range(0,len(snaparr[i])):
			for k in range(0,len(T_ea)):
#
	# Define which EA snapshots occur during outburst phase
#
				outburst = 0
				for a in range(0, n_accr):
					if ((T_ea[k] > 1) and (T_time[k] > t_s[a]) and (T_time[k] < t_e[a])):
						outburst = 1
#
				if ((snaparr[i][j] == T_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 0) and (outburst == 0)):
					plt.scatter(T_time[k], abs(T_exp[k]), s=80, \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
					ax1.set_yticks(ytck)
					n_out = n_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 1) and (outburst == 0)):
					plt.scatter(T_time[k], abs(T_exp[k]), s=80, marker="s", \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
					ax1.set_yticks(ytck)
					n_out = n_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (T_ea[k] == ea_run_n[i]) and (exp_err[k] == 0) and (outburst == 1)):
					plt.scatter(T_time[k], abs(T_exp[k]), s=80, \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
					ax1.set_yticks(ytck)
					y_out = y_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (T_ea[k] == ea_run_n[i]) and (exp_err[k] == 1) and (outburst == 1)):
					plt.scatter(T_time[k], abs(T_exp[k]), s=80, marker="s", \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
					ax1.set_yticks(ytck)
					y_out = y_out + 1
				plt.ylim(0, 2.5) ; plt.xlim(78., 100.) ; plt.ylabel('q')
			plt.legend(loc='upper right', fontsize = 8, scatterpoints = 1)
	plt.xlabel('Time (kyr)' )
#
	plt.savefig(plotdir+'T_exp.pdf') ; plt.clf()

