#
# pdisc_r.py
#
# Python program to read in pdisc (radially varied, temporally static) data structures. 
# Uses information of file name to generate plots analogous to Stamatellos, 
# Whitworth & Hubber (2012) work.
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
	# !!! DO TIME CHECK PRIOR TO ENTERING VALUES OF [tag]_s_[time] INDICES !!!
#
powerfit_check = "FALSE"		# Choose whether ("TRUE") or not ("FALSE") to fitted power indices to Sig and T profiles
#
col_arr = ["b", "g", "r", "c", "m", "k"]
#
#
def func(x, a, b):			# Linear function used to fit power indices to log-log Sig/T vs. radius distributions
  return a*x + b
#
smooth = 3				# Set level of smoothing that SD and  distributions undergo (1 = original data)
loglin = "FALSE"				# Choose either log-linear ("TRUE") or log-log ("FALSE") distributions of SD and T
#
r_start = 10				# Radial index to start plots (must be > 1 for loglin == "FALSE")
#
d = {}						# Radial value (AU) to start of p and q fits for each snaparr (see main.py for formatting)
r_fitS_0 = 'r_fitS_0' ; r_fitS_1 = 'r_fitS_1'	# using dictionaries
r_fitS_3 = 'r_fitS_3' ; r_fitS_4 = 'r_fitS_4' ; r_fitS_5 = 'r_fitS_5'	
d[r_fitS_0] = [25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25]
d[r_fitS_1] = [25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25]
d[r_fitS_3] = [[20,20,20,20],[20,20,20,20],[20,20,20,20]]
d[r_fitS_4] = [[15,15,20,15],[20,20,20,15],[20,20,20,20]]
d[r_fitS_5] = [[15,15,15,15],[15,15,15,15],[15,15,15,15]]
#
r_fitE_0 = 'r_fitE_0' ; r_fitE_1 = 'r_fitE_1'	# As above, for end of fits for each snaparr.
r_fitE_3 = 'r_fitE_3' ; r_fitE_4 = 'r_fitE_4' ; r_fitE_5 = 'r_fitE_5'	
d[r_fitE_0] = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50]
d[r_fitE_1] = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50]
d[r_fitE_3] = [[50,50,50,50],[50,50,50,50],[50,50,50,50]]
d[r_fitE_4] = [[40,35,40,35],[50,50,50,40],[60,60,60,60]]
d[r_fitE_5] = [[40,40,40,40],[60,60,60,40],[60,70,70,70]]
#
exp_err_0 = 'exp_err_0' ; exp_err_1 = 'exp_err_1'	# Error tag on fitting to SD profile, to indicate "bumps" in fitting routine
exp_err_3 = 'exp_err_3' ; exp_err_4 = 'exp_err_4' ; exp_err_5 = 'exp_err_5'	
d[exp_err_0] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
d[exp_err_1] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
d[exp_err_3] = [[0, 0, 0, 0],[0, 0, 0, 1],[1, 0, 0, 0]]
d[exp_err_4] = [[0, 0, 1, 0],[0, 1, 1, 1],[0, 0, 0, 0]]
d[exp_err_5] = [[0, 1, 0, 0],[0, 0, 0, 1],[1, 1, 1, 1]]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
def read(dat_dir, plotdir, ea_run, snaparr, snapcore, EA_timeref, EA_lenref, pmass, \
   pradius, hasharr_app, n_accr, r_limit, r_d_kep, r_d_sigALMA, timearr, v_K, inclin):
#
	print("pdisc files being read")
#
	# Define number of files to be read dependent on snaparr dimensions
#
	if (snaparr.ndim == 1):
		file_n = len(snaparr)
	elif (snaparr.ndim == 2):
		file_n = len(snaparr)*len(snaparr[0])
#
	# Arrays to fill
#
	time = [[] for i in range(file_n)] ; r = [[] for i in range(file_n)]
	Q = [[] for i in range(file_n)] ; T = [[] for i in range(file_n)]
	sig = [[] for i in range(file_n)] ; rInM = [[] for i in range(file_n)]
	M_s = [[] for i in range(file_n)] ; M_d = [[] for i in range(file_n)]
	v_r = [[] for i in range(file_n)] ; v_the = [[] for i in range(file_n)]
	v_z = [[] for i in range(file_n)] ; vkep = [[] for i in range(file_n)]
	v_mod = [[] for i in range(file_n)]
#
	# Create string array as pointer to pdisc files based on snaparr
#
	file_list = []
	if (snaparr.ndim == 1):
		snaparr_tmp = np.array([0]*len(snaparr))
		snapcore_tmp = np.array([0]*len(snapcore))
		fcount = 0
		for a in range(0, len(snaparr)):
			if (snaparr[a] < (1000-70)):
				file_list.append(dat_dir+'pdisc/DE05.du.00'+ \
				   str(snaparr[a]+69)+'.pdisc.1')
			elif (snaparr[a] > (1000-70)):
				file_list.append(dat_dir+'pdisc/DE05.du.0'+ \
				   str(snaparr[a]+69)+'.pdisc.1')
			for b in range(0, len(snapcore)):
				if (snaparr[a] == snapcore[b]):
					snapcore_tmp[b] = fcount
			snaparr_tmp[a] = fcount
			fcount = fcount + 1
#
	elif (snaparr.ndim == 2):
		snaparr_tmp = np.array([[0]*len(snaparr[0])]*len(snaparr)) 
		fcount = 0
		for a in range(0, len(snaparr)):
			for b in range(0, len(snaparr[0])):
				if (snaparr[a][b] < (1000-70)):
					file_list.append(dat_dir+'pdisc/DE05.du.00'+ \
					   str(snaparr[a][b]+69)+'.pdisc.1')
				elif (snaparr[a][b] > (1000-70)):
					file_list.append(dat_dir+'pdisc/DE05.du.0'+ \
					   str(snaparr[a][b]+69)+'.pdisc.1')
				snaparr_tmp[a][b] = fcount
				fcount = fcount + 1
#
	# Loop over files in pdisc directory to read all temporally evolved values
#
	for i in range(0, file_n):
#
	# Now loop over lines in each file and extract variables of interest
	# then convert to numpy format for future data extraction for plotting
#
		f = open(file_list[i], 'r')
#
	# Firstly, read and ignore file header
#
		header = f.readline()
#
		for line in f:
			line = line.strip() ; columns = line.split()
			time[i].append(float(columns[0])/1000.) ; r[i].append(float(columns[1]))
			Q[i].append(float(columns[2])) ; T[i].append(float(columns[3]))
			sig[i].append(float(columns[4])) ; rInM[i].append(float(columns[7]))
			M_s[i].append(float(columns[8])) ; M_d[i].append(float(columns[9]))
			v_r[i].append(float(columns[10])) ; v_the[i].append(float(columns[11]))
			v_z[i].append(float(columns[12])) ; vkep[i].append(float(columns[13]))
			v_mod[i].append(math.sqrt(float(columns[10])**2 + float(columns[11])**2 \
			   + float(columns[12])**2))
		f.close()
#
	# Set run specific fit start locations as per dictionary definitions
#
		r_fit_S = d['r_fitS_'+str(ea_run)] ; r_fit_S = np.array(r_fit_S)
		r_fit_E = d['r_fitE_'+str(ea_run)] ; r_fit_E = np.array(r_fit_E)
		exp_err = d['exp_err_'+str(ea_run)]
#
	time = np.array(time) ; r = np.array(r) ; Q = np.array(Q) ; T = np.array(T) ; sig = np.array(sig)
	rInM = np.array(rInM) ; M_s = np.array(M_s) ; M_d = np.array(M_d) ; v_r = np.array(v_r)
	v_the = np.array(v_the) ; v_z = np.array(v_z) ; vkep = np.array(vkep) ; v_mod = np.array(v_mod)
#
#
#
	# Plot the disc parameters before, during and after the accretion for defined radial value
#
	print("pdisc plots now being generated")
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# - - - NON-EA RUNS
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
# Begin plotting comparison of simulation time references
#
	if (snaparr_tmp.ndim == 1):
#
		fig = plt.figure(1)
#
		fig.subplots_adjust(hspace=.3, wspace = 0.4)
		fig.suptitle('Radially varied disc parameters for R < '+str(r_limit)+' AU', \
		   fontsize=10)
#
		ax1 = plt.subplot(221)
		for i in range(0, len(snaparr_tmp)):
			for j in range(0, len(snapcore_tmp)):
				if (snaparr_tmp[i] == snapcore_tmp[j]):
					line1 = plt.plot(r[snapcore_tmp[j],0:r_limit], \
					Q[snapcore_tmp[j],0:r_limit], color = col_arr[j])
		plt.xlabel("Disc radius (AU)", fontsize = 12) 
		plt.ylabel("Toomre Q parameter", fontsize = 12, labelpad=0.5)
		plt.ylim(0, 4)
		plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
#
		ax2 = plt.subplot(222)
		for i in range(0, len(snaparr_tmp)):
			for j in range(0, len(snapcore_tmp)):
				if (snaparr_tmp[i] == snapcore_tmp[j]):
					line2 = plt.plot(r[snapcore_tmp[j],0:r_limit], \
					   np.log10(T[snapcore_tmp[j],0:r_limit]), \
					   label = "t = "+str(round(timearr[i], 2))+" kyr ", \
					   color = col_arr[j])
		plt.xlabel("Disc radius (AU)", fontsize = 12) 
		plt.ylabel("log Temperature (K)", fontsize = 12, labelpad=0.5)
		plt.ylim(ax2.get_ylim()[0], ax2.get_ylim()[1])
 		plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
		plt.legend(loc='upper right', fontsize = 10)
#
		ax3 = plt.subplot(223)
		for i in range(0, len(snaparr_tmp)):
			for j in range(0, len(snapcore_tmp)):
				if (snaparr_tmp[i] == snapcore_tmp[j]):
					line3 = plt.plot(r[snapcore_tmp[j],0:r_limit], \
					   rInM[snapcore_tmp[j],0:r_limit], color = col_arr[j])
		plt.xlabel("Disc radius (AU)", fontsize = 12) 
		plt.ylabel("Disc mass "+(r'(M$_{\odot}$)'), fontsize = 12, labelpad=0.5)
		plt.ylim(ax3.get_ylim()[0], ax3.get_ylim()[1])
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
#
		ax4 = plt.subplot(224)
		for i in range(0, len(snaparr_tmp)):
			for j in range(0, len(snapcore_tmp)):
				if (snaparr_tmp[i] == snapcore_tmp[j]):
					line4 = plt.plot(r[snapcore_tmp[j],0:r_limit], \
					   v_r[snapcore_tmp[j],0:r_limit], color = col_arr[j])
		plt.xlabel("Disc radius (AU)", fontsize = 12) 
		plt.ylabel((r'v$_{r}$')+' (km'+(r's$^{-1}$')+')', fontsize = 12, labelpad=0.5)
		plt.ylim(ax4.get_ylim()[0], ax4.get_ylim()[1])
		plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
#
		plt.savefig(str(plotdir)+'Matrix_'+str(r_limit)+'AU.pdf') ; plt.clf()
#
	# Plot individual r vs. Q distribution
#
		fig = plt.figure(1)
		ax1 = plt.subplot(111)
		for i in range(0, len(snaparr_tmp)):
			for j in range(0, len(snapcore_tmp)):
				if (snaparr_tmp[i] == snapcore_tmp[j]):
					line1 = plt.plot(r[snapcore_tmp[j],0:r_limit], \
					   Q[snapcore_tmp[j],0:r_limit],
					   label = "t = "+str(round(timearr[i], 2))+" kyr ", \
					   color = col_arr[j])
		plt.xlabel("Disc radius (AU)", fontsize = 12) 
		plt.ylabel("Toomre Q parameter", fontsize = 12, labelpad=0.5)
		ax1.set_ylim(0, 4)
		plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
		plt.legend(loc='upper right', fontsize = 10)
		plt.savefig(str(plotdir)+'Q_r_'+str(r_limit)+'AU.pdf') ; plt.clf()
#
# Plot surface density vs. radius for non-ea runs with fit to log-log data to find power index, p
#
	# Set y_fitted array sizes dependent on disc radial extents using proxy of ALMA sig. criterion
# 
		coeffs = [0]*len(snaparr_tmp) ; matcov = [0]*len(snaparr_tmp)
		y_fitted = [[] for i in range(len(snaparr_tmp))]
#
		for i in range(0, len(snaparr_tmp)):
			for j in range(0, len(sig[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]] ) ):
				y_fitted[i].append(0)
#
		plt.figure(1)
		ax1 = plt.subplot(111)
#
		for i in range(0, len(snaparr_tmp)):
#
			coeffs[i], matcov[i] = curve_fit(func, np.log10(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), \
			   np.log10(sig[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), [1, 1])
#
	# Write exponents for SD
#
			if( (v_K == "90") and (inclin == "0")):
				f = open(dat_dir+'../../SD_exps.dat','a')
				f.write(str(ea_run)+' '+str(snaparr[i])+' '+str(timearr[i])+' ' \
				   +str( round(coeffs[i][0], 2) )+' '+str(exp_err[i])+'\n' )
				f.close()
#
			for j in range(0, len(snapcore_tmp)):
				if (snaparr_tmp[i] == snapcore_tmp[j]):
#
					y_fitted[i] = func(np.log10(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), \
					    coeffs[i][0], coeffs[i][1])
#
					rnew = np.arange(r_start,r_limit,r_limit/(r_limit/smooth))
#
					tck = interpolate.splrep(r[snaparr_tmp[i],r_start:r_limit], \
					   sig[snaparr_tmp[i],r_start:r_limit])
					signew = interpolate.splev(rnew, tck, der = 0)
#
					if (loglin == "TRUE"):
						line1 = plt.plot(rnew, signew, color = col_arr[j])
						line2 = plt.plot(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]], \
						   np.power(10.,y_fitted[i]), color = col_arr[j], linewidth = 2, \
						   linestyle = 'dashed', label = str(round(coeffs[i][0], 1)) )
						ax1.set_yscale('log')
					elif (loglin == "FALSE"):
						line1 = plt.plot(rnew, signew, color = col_arr[j])
						line2 = plt.plot(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]], \
						   np.power(10.,y_fitted[i]), color = col_arr[j], linewidth = 2, \
						   linestyle = 'dashed', label = str(round(coeffs[i][0], 1)) )
						ax1.set_xscale('log') ; ax1.set_yscale('log')
#
		plt.xlabel("Disc Radius (AU)", fontsize = 12) 
		plt.ylabel("log Surface density", fontsize = 12, labelpad=0.5)
		ax1.set_ylim(0, 1550) ; ax1.set_xlim(0, 155)
		plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
		plt.legend(loc='upper right', title = "p values", fontsize=10)
#
		plt.savefig(str(plotdir)+'SD_r_'+str(r_limit)+'AU.pdf') ; plt.clf()	
#
	# Print power law fits to surface density profiles
#
		if (powerfit_check == "TRUE"):
			print "\n Surface density power law fits \n"
			for i in range(0, len(snaparr_tmp)):
				print "Power index for the "+snaparr_tmp[i]+" snapshot of " \
				   +str(ea_run)+" is: ",  \
				   str( round( coeffs[i][0] , 2 ) )+" \n"
#
# Plot temperature vs. radius for non-ea runs with fit to log-log data to find power index, q
#
	# Set y_fitted array sizes dependent on disc radial extents using proxy of ALMA sig. criterion
# 
		coeffs = [0]*len(snaparr_tmp) ; matcov = [0]*len(snaparr_tmp)
		y_fitted = [[] for i in range(len(snaparr_tmp))]
#
		for i in range(0, len(snaparr_tmp)):
			for j in range(0, len(T[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]] ) ):
				y_fitted[i].append(0)
#
		plt.figure(1)
		ax1 = plt.subplot(111)
		for i in range(0, len(snaparr_tmp)):
#
			coeffs[i], matcov[i] = curve_fit(func, np.log10(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), \
			   np.log10(T[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), [1, 1])
			y_fitted[i] = func(np.log10(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), \
			   coeffs[i][0], coeffs[i][1])
#
	# Write exponents for T
#
			if( (v_K == "90") and (inclin == "0")):
				f = open(dat_dir+'../../T_exps.dat','a')
				f.write(str(ea_run)+' '+str(snaparr[i])+' '+str(timearr[i])+ \
				   ' '+str( round(coeffs[i][0], 2) )+'\n' )
				f.close()
#
			for j in range(0, len(snapcore_tmp)):
				if (snaparr_tmp[i] == snapcore_tmp[j]):
#
					rnew = np.arange(r_start,r_limit,r_limit/(r_limit/smooth))
#
					tck = interpolate.splrep(r[snaparr_tmp[i],r_start:r_limit], \
					   T[snaparr_tmp[i],r_start:r_limit])
					Tnew = interpolate.splev(rnew, tck, der = 0)
#
					if (loglin == "TRUE"):
						line1 = plt.plot(rnew, Tnew, color = col_arr[j])
						line2 = plt.plot(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]], \
						   np.power(10.,y_fitted[i]), color = col_arr[j], linewidth = 2, \
						   linestyle = 'dashed', label = str(round(coeffs[i][0], 1)) )
						ax1.set_yscale('log')
#
					elif (loglin == "FALSE"):
						line1 = plt.plot(rnew, Tnew, color = col_arr[j])
						line2 = plt.plot(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]], \
						   np.power(10,y_fitted[i]), color = col_arr[j], linewidth = 2, \
						   linestyle = 'dashed', label = str(round(coeffs[i][0], 1)) )
						ax1.set_xscale('log') ; ax1.set_yscale('log')
#
		plt.xlabel("Disc Radius (AU)", fontsize = 12) 
		plt.ylabel("log Temperature (K)", fontsize = 12, labelpad=0.5)
		ax1.set_ylim(0, 550) ; ax1.set_xlim(0, 155)
		plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
		plt.legend(loc='upper right', title = "q values", fontsize=10)
#
		plt.savefig(str(plotdir)+'T_r_'+str(r_limit)+'AU.pdf') ; plt.clf()	
#
	# Print power law fits to Temperature profiles
#
		if (powerfit_check == "TRUE"):
			print "\n Temperature power law fits \n"
			for i in range(0, len(snaparr_tmp)):
				print "Power index for the "+snaparr_tmp[i]+" snapshot of "+ \
				   str(ea_run)+" is: ",  \
				   str(round( coeffs[i][0], 2 ) )+" \n"
#
# Plot azimuthal velocity vs. radius for singular accretion event. Compare data of before/during/
# after event with analytic calculation of Keplerian velocity for snapshot
#
		line1 = [ []*len(v_the[0,0:r_limit]) ]*len(snapcore_tmp)
		line2 = [ []*len(v_the[0,0:r_limit]) ]*len(snapcore_tmp)
		plt.figure(1)
		ax1 = plt.subplot(111)
		plt.title('Azimuthal velocity vs. radius')
		for i in range(0, len(snaparr_tmp)):
			for j in range(0, len(snapcore_tmp)):
				if (snaparr_tmp[i] == snapcore_tmp[j]):
					line1[j] = plt.plot(r[snaparr_tmp[i],0:r_limit], \
					   abs(v_the[snaparr_tmp[i],0:r_limit]), color = col_arr[j])
					line2[j] = plt.plot(r[snaparr_tmp[i],0:r_limit], \
					   vkep[snaparr_tmp[i],0:r_limit], \
					   label = 'Keplerian velocity ('+str(snaparr_tmp[i])+' analytical)', \
					   linewidth = 2, linestyle = 'dashed', color=col_arr[j])
					plt.axvline(x=r_d_kep[snaparr_tmp[i]], linewidth = 2, \
					   linestyle = 'dotted', color=col_arr[j])
#
	# If a planet is present, overplot onto x-axis with hashed line
#
		if (len(pmass) > 0):
			for a in range(1, len(pmass)):
				for idash in range(1, len(pmass[0])):
					if (pmass[a][idash] != 0.):
						ax1.plot(pradius[a][idash], \
						   vkep[snaparr_tmp[idash],int(pradius[a][idash])], \
						    col_arr[idash]+'o', markersize=10.0)
		plt.xlabel("Disc Radius (AU)", fontsize = 12)
		plt.ylabel(''+(r'v$_{\phi}$')+' (km'+(r's$^{-1}$')+')', fontsize = 12, labelpad=0.5)
		plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1]/2.) ; plt.xlim(0, r_limit)
		plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
		plt.legend(loc='upper right', fontsize = 10)
#
		plt.savefig(str(plotdir)+'Azi_r_'+str(r_limit)+'AU.pdf') ; plt.clf()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# - - - EA RUNS [ VARYING OUTBURST DURATION ]
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	if (snaparr_tmp.ndim == 2):
#
	# Invoke a switch statement to loop over pdisc analysis of runs with ea 
	# to focus on (a) single accretion event and (b) focus on snapshot relative to accretion
	# state for different events
#	
		for lendur in range(0, 2):
			title_point = [] ; leg_point = []
#
			if (lendur == 0):
				for i in range(0, len(snaparr_tmp)):
					title_point.append(EA_lenref[i])
				for i in range(0, len(snaparr_tmp[0])):
					leg_point.append(EA_timeref[i])
			elif (lendur == 1):
				for i in range(0, len(snaparr_tmp[0])):
					title_point.append(EA_timeref[i])
				for i in range(0, len(snaparr_tmp)):
					leg_point.append(EA_lenref[i])
#
	# Rotate snaparr, timearr and planet mass/radius data dependent on switch
#
				snaparr_tmp = np.rot90(snaparr_tmp, 1) ; snaparr_tmp = snaparr_tmp[::-1]
				timearr = np.rot90(timearr, 1) ; timearr = timearr[::-1]
				r_fit_S = np.rot90(r_fit_S, 1) ; r_fit_S = r_fit_S[::-1]
				r_fit_E = np.rot90(r_fit_E, 1) ; r_fit_E = r_fit_E[::-1]
				for a in range(1, len(pmass)):
					pmass[a] = np.rot90(pmass[a], 1) ; pmass[a] = pmass[a][::-1]
					pradius[a] = np.rot90(pradius[a], 1) ; pradius[a] = pradius[a][::-1]
#
			for i in range(0,len(title_point)):
#
				fig = plt.figure(1)
#
				fig.subplots_adjust(hspace=.3, wspace = 0.4)
				fig.suptitle('Disc parameters '+title_point[i]+' for R < '+str(r_limit)+ \
				   ' AU', fontsize=10)
#
				ax1 = plt.subplot(221)
				for j in range(0, len(snaparr_tmp[0])):
					line1 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
					   Q[snaparr_tmp[i][j],0:r_limit], color = col_arr[j])
				plt.xlabel("Disc radius (AU)", fontsize = 12) 
				plt.ylabel("Toomre Q parameter", fontsize = 12, labelpad=0.5)
				plt.ylim(0, 4)
				plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
#
				ax2 = plt.subplot(222)
				for j in range(0, len(snaparr_tmp[0])):
					line2 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
					   np.log10(T[snaparr_tmp[i][j],0:r_limit]), \
					   label = "t = "+str(round(timearr[i][j], 2))+" kyr "+ \
					   EA_timeref[j], color = col_arr[j])
				plt.xlabel("Disc radius (AU)", fontsize = 12) 
				plt.ylabel("log Temperature (K)", fontsize = 12, labelpad=0.5)
				plt.ylim(ax2.get_ylim()[0], ax2.get_ylim()[1])
 				plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
				plt.legend(loc='upper right', fontsize = 10)
				ax1.set_yscale('log')
#
				ax3 = plt.subplot(223)
				for j in range(0, len(snaparr_tmp[0])):
					line3 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
					   rInM[snaparr_tmp[i][j],0:r_limit], color = col_arr[j])
				plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
				plt.xlabel("Disc radius (AU)", fontsize = 12) 
				plt.ylabel("Disc mass "+(r'(M$_{\odot}$)'), fontsize = 12, labelpad=0.5)
				plt.ylim(ax3.get_ylim()[0], ax3.get_ylim()[1])
				plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
#
				ax4 = plt.subplot(224)
				for j in range(0, len(snaparr_tmp[0])):
					line4 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
					   v_r[snaparr_tmp[i][j],0:r_limit], color = col_arr[j])
				plt.xlabel("Disc radius (AU)", fontsize = 12) 
				plt.ylabel((r'v$_{r}$')+' (km'+(r's$^{-1}$')+')',fontsize = 12,labelpad=0.5)
				plt.ylim(ax4.get_ylim()[0], ax4.get_ylim()[1])
				plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
#
				plt.savefig(str(plotdir)+'Matrix_'+str(r_limit)+'AU_'+title_point[i]+'.pdf')
				plt.clf()
#
	# Plot individual r vs. Q distribution
#
				fig = plt.figure(1)
				ax1 = plt.subplot(111)
				for j in range(0, len(snaparr_tmp[0])):
					line1 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
					   Q[snaparr_tmp[i][j],0:r_limit], color = col_arr[j])
				plt.xlabel("Disc radius (AU)", fontsize = 12) 
				plt.ylabel("Toomre Q parameter", fontsize = 12, labelpad=0.5)
				ax1.set_ylim(0, 4)
				plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
				plt.savefig(str(plotdir)+'Q_r_'+str(r_limit)+'AU_'+title_point[i]+'.pdf')
				plt.clf()
#
# Plot surface density vs. radius for ea runs with fit to log-log data to find power index, p
#
				coeffs = [0]*len(snaparr_tmp[0]) ; matcov = [0]*len(snaparr_tmp[0])
				y_fitted = [[] for j in range(len(snaparr_tmp[0]))]
#
	# Set y_fitted array sizes dependent on disc radial extents using proxy of ALMA sig. criterion
# 
				for j in range(0, len(snaparr_tmp[0])):
#
					for k in range(0, len(sig[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]) ):
						y_fitted[j].append(0)
#
				plt.figure(1)
				ax1 = plt.subplot(111)
#
				for j in range(0, len(snaparr_tmp[0])):
#
					coeffs[j], matcov[j] = curve_fit(func, np.log10(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), \
					   np.log10(sig[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), [1, 1])
					y_fitted[j] = func(np.log10(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), \
					   coeffs[j][0], coeffs[j][1])
#
	# Write exponents for SD
#
					if( (v_K == "90") and (inclin == "0") and (lendur == 0) ):
						f = open(dat_dir+'../../SD_exps.dat','a')
						f.write(str(ea_run)+' '+str(snaparr[i][j])+' '+str(timearr[i][j])+ \
						   ' '+str( round(coeffs[j][0], 2) )+' '+str(exp_err[i][j])+'\n' )
						f.close()
#
					rnew = np.arange(r_start,r_limit,r_limit/(r_limit/smooth))
#
					tck = interpolate.splrep(r[snaparr_tmp[i][j],r_start:r_limit], \
					   sig[snaparr_tmp[i][j],r_start:r_limit])
					signew = interpolate.splev(rnew, tck, der = 0)
#
					if (loglin == "TRUE"):
						line1 = plt.plot(rnew, signew, color = col_arr[j])
						line2 = plt.plot(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]], \
						   np.power(10.,y_fitted[j]), color = col_arr[j], linewidth = 2, \
						   linestyle = 'dashed', \
						   label = str(round(coeffs[j][0], 1)) )
						ax1.set_yscale('log')
#
					elif (loglin == "FALSE"):
						line1 = plt.plot(rnew, signew, color = col_arr[j])
						line2 = plt.plot(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]], \
						   np.power(10.,y_fitted[j]), color = col_arr[j], linewidth = 2, \
						   linestyle = 'dashed', \
						   label = str(round(coeffs[j][0], 1)) )
						ax1.set_xscale('log') ; ax1.set_yscale('log')
#
				plt.xlabel("Disc Radius (AU)", fontsize = 12) 
				plt.ylabel("log Surface Density", fontsize = 12, labelpad=0.5)
				ax1.set_ylim(0, 1250) ; ax1.set_xlim(0, 155)
				plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
				plt.legend(loc='upper right', title = "p values", \
				   fontsize=10)
#	
				plt.savefig(str(plotdir)+'SD_r_'+str(r_limit)+'AU_'+str(title_point[i])+'.pdf')
				plt.clf()	
#
	# Print power law fits to surface density profiles
#
				if (powerfit_check == "TRUE"):
					print "\n Surface density power law fits \n"
					for j in range(0, len(snaparr_tmp)):
						print "Power index for the "+str(snaparr_tmp[i][j]) \
						   +" snapshot of "+str(ea_run)+" is: ",  \
						   str( round( coeffs[j][0] , 2 ) )+" \n"
#
# Plot temperature vs. radius for ea runs with fit to log-log data to find power index, q
#
				coeffs = [0]*len(snaparr_tmp[0]) ; matcov = [0]*len(snaparr_tmp[0])
				y_fitted = [[] for j in range(len(snaparr_tmp[0]))]
#
	# Set y_fitted array sizes dependent on disc radial extents using proxy of ALMA sig. criterion
# 
				for j in range(0, len(snaparr_tmp[0])):
					for k in range(0, len(T[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]) ):
						y_fitted[j].append(0)
#
				plt.figure(1)
				ax1 = plt.subplot(111)
#
				for j in range(0, len(snaparr_tmp[0])):
#
					coeffs[j], matcov[j] = curve_fit(func, np.log10(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), \
					   np.log10(T[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), [1, 1])
					y_fitted[j] = func(np.log10(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), \
					   coeffs[j][0], coeffs[j][1])
#
	# Write exponents for T
#
					if( (v_K == "90") and (inclin == "0") and (lendur == 0) ):
						f = open(dat_dir+'../../T_exps.dat','a')
						f.write(str(ea_run)+' '+str(snaparr[i][j])+' '+str(timearr[i][j])+ \
						   ' '+str( round(coeffs[j][0], 2) )+'\n' )
						f.close()
#
					rnew = np.arange(r_start,r_limit,r_limit/(r_limit/smooth))
#
					tck = interpolate.splrep(r[snaparr_tmp[i][j],r_start:r_limit], \
					   T[snaparr_tmp[i][j],r_start:r_limit])
					Tnew = interpolate.splev(rnew, tck, der = 0)
#
					if (loglin == "TRUE"):
						line1 = plt.plot(rnew, Tnew, color = col_arr[j])
						line2 = plt.plot(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]], \
						   np.power(10.,y_fitted[j]), color = col_arr[j], linewidth = 2, \
						   linestyle = 'dashed', \
						   label = str(round(coeffs[j][0], 1)) )
						ax1.set_yscale('log')
#
					if (loglin == "FALSE"):
						line1 = plt.plot(rnew, Tnew, color = col_arr[j])
						line2 = plt.plot(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]], \
						   np.power(10.,y_fitted[j]), color = col_arr[j], linewidth = 2, \
						   linestyle = 'dashed', \
						   label = str(round(coeffs[j][0], 1)) )
						ax1.set_xscale('log') ; ax1.set_yscale('log')
#
				plt.xlabel("Disc Radius (AU)", fontsize = 12) 
				plt.ylabel("log Temperature (K)", fontsize = 12, labelpad=0.5)
				ax1.set_ylim(0, 350) ; ax1.set_xlim(0, 155)
				plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
				plt.legend(loc='upper right', title = "q values", \
				   fontsize=10)
#
				plt.savefig(str(plotdir)+'T_r_'+str(r_limit)+'AU_'+str(title_point[i])+'.pdf')
				plt.clf()	
#
	# Print power law fits to surface density profiles
#
				if (powerfit_check == "TRUE"):
					print "\n Temperature power law fits \n"
					for j in range(0, len(snaparr_tmp[0])):
						print "Power index for the "+str(snaparr_tmp[i][j]) \
						   + " snapshot of "+str(ea_run)+" is: ",  \
						   str( round( coeffs[j][0] , 2 ) ) +" \n"
#
# Compare velocity data with analytic calculation of Keplerian velocity
#
				plt.figure(1)
				ax1 = plt.subplot(111)
				plt.title('Azimuthal velocity vs. radius for '+str(title_point[i])+' run')
#
				for j in range(0, len(leg_point)):

					line1 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
					   abs(v_the[snaparr_tmp[i][j],0:r_limit]), color=col_arr[j])
					line2 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], 
					   vkep[snaparr_tmp[i][j],0:r_limit], \
					   label = 'Keplerian velocity ('+str(leg_point[j])+' analytical)', \
					   linewidth = 2, linestyle = 'dashed', color=col_arr[j])
#
	# If a planet is present, overplot onto computed Keplerian profile with filled circle
#
					if (len(pmass) > 1):
						for a in range(1, len(pmass)):
							if (pmass[a][i][j] != 0.):
								ax1.plot(pradius[a][i][j],
								   vkep[snaparr_tmp[i][j],int(pradius[a][i][j])], \
								   col_arr[j]+'o', markersize=10.0)
#
					if (r_d_kep != 0):
						plt.axvline(x=r_d_kep[snaparr_tmp[i][j]], linewidth = 2, \
						   linestyle = 'dotted', color = col_arr[j])
#
				plt.xlabel("Disc radius (AU)", fontsize = 12)
				plt.ylabel(''+(r'v$_{\phi}$')+' (km'+(r's$^{-1}$')+')', \
				   fontsize = 12, labelpad=0.5)
				plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1]/2.) ; plt.xlim(0, r_limit)
				plt.yticks(fontsize = 12) ; plt.xticks(fontsize = 12)
				plt.legend(loc='upper right', fontsize = 10)
#
				plt.savefig(str(plotdir)+'Azi_r_'+str(r_limit)+'AU_'+str(title_point[i])+'.pdf')
				plt.clf()
#
	# Derotate snaparr_tmp, timearr and planet mass/radius arrays for further analyses
#
#
		snaparr_tmp = snaparr_tmp[::-1] ; snaparr_tmp = np.rot90(snaparr_tmp, -1)
		timearr = timearr[::-1] ; timearr = np.rot90(timearr, -1)
		r_fit_S = r_fit_S[::-1] ; r_fit_S = np.rot90(r_fit_S, -1) 
		r_fit_E = r_fit_E[::-1] ; r_fit_E = np.rot90(r_fit_E, -1) 
		for a in range(1, len(pmass)):
			pmass[a] = np.array(pmass[a][::-1]) ; pmass[a] = np.rot90(pmass[a], -1)
			pradius[a] = np.array(pradius[a][::-1]) ; pradius[a] = np.rot90(pradius[a], -1)
#
#
	return r, vkep, col_arr
