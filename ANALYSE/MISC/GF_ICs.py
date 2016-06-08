#
# GI_ICs.py
#
# Python program to read in GI initial conditions from selected pdisc data structures 
# from SEREN analysedisc executable. Uses information of file name to generate plots 
# analogous to Stamatellos, Whitworth & Hubber (2012) work
#
# Author: Benjamin MacFarlane
# Date: 07/03/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit
from scipy import interpolate
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# !!! DO TIME CHECK PRIOR TO ENTERING VALUES OF [tag]_s_[time] INDICES !!!
#
#
arch_dir = "/home/ben/Documents/WORK_PLANETS/PROJECTS/FUORS/DATA/3/GI_ICs/"	# Location of ea/ directory with data set
plotdir = arch_dir+"../../../PLOTS/GI_ICs/"
#
ea_run = 1
#
r_start = 1
r_limit = 100
#
spline = "TRUE"
spline_smooth = 3
#
snaparr = [606,1356]
#
powerfit_check = "TRUE"
#
col_arr = ["b", "g"]
#
def func(x, a, b):
  return a*x + b
#
r_fit = 5				# In AU
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
os.system('mkdir '+str(plotdir))
#
	# Define arrays to fill
#
file_n = len(snaparr)
#
time = [[] for i in range(file_n)] ; r = [[] for i in range(file_n)]
Q = [[] for i in range(file_n)] ; T = [[] for i in range(file_n)]
sig = [[] for i in range(file_n)] ; Omeg = [[] for i in range(file_n)]
DiscM = [[] for i in range(file_n)] ; rInM = [[] for i in range(file_n)]
M_s = [[] for i in range(file_n)] ; M_d = [[] for i in range(file_n)]
v_r = [[] for i in range(file_n)] ; v_the = [[] for i in range(file_n)]
v_z = [[] for i in range(file_n)] ; vkep = [[] for i in range(file_n)]
h = [[] for i in range(file_n)] ; h_mid = [[] for i in range(file_n)]
rpart = [[] for i in range(file_n)] ; v_mod = [[] for i in range(file_n)]
#
	# Must remake snaparr array to point to refined pdisc files read in
#
file_list = []
snaparr_tmp = np.array([0]*len(snaparr))
fcount = 0
for a in range(0, len(snaparr)):
	if (snaparr[a] < (1000-70)):
		file_list.append(arch_dir+'DE05.du.00'+ \
		   str(snaparr[a]+69)+'.pdisc.1')
	elif (snaparr[a] > (1000-70)):
		file_list.append(arch_dir+'DE05.du.0'+ \
		   str(snaparr[a]+69)+'.pdisc.1')
	snaparr_tmp[a] = fcount
	fcount = fcount + 1
#
#
	# Loop over files in pdisc directory to read all temporally evolved values
#
for i in range(0,file_n):
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
		sig[i].append(float(columns[4])) ; Omeg[i].append(float(columns[5]))
		DiscM[i].append(float(columns[6])) ; rInM[i].append(float(columns[7]))
		M_s[i].append(float(columns[8])) ; M_d[i].append(float(columns[9]))
		v_r[i].append(float(columns[10])) ; v_the[i].append(float(columns[11]))
		v_z[i].append(float(columns[12])) ; vkep[i].append(float(columns[13]))
		h[i].append(float(columns[14])) ; h_mid[i].append(float(columns[15]))
		rpart[i].append(float(columns[16]))
		v_mod[i].append(math.sqrt(float(columns[10])**2 + float(columns[11])**2 + \
		   float(columns[12])**2))
	f.close()
#
time = np.array(time) ; r = np.array(r) ; Q = np.array(Q) ; T = np.array(T)
sig = np.array(sig) ; Omeg = np.array(Omeg) ; DiscM = np.array(DiscM)
rInM = np.array(rInM) ; M_s = np.array(M_s) ; M_d = np.array(M_d)
v_r = np.array(v_r) ; v_the = np.array(v_the) ; v_z = np.array(v_z)
vkep = np.array(vkep) ; h = np.array(h) ; h_mid = np.array(h_mid)
rpart = np.array(rpart) ; v_mod = np.array(v_mod)
#
# Plot surface density vs. radius for singular accretion event
# First, fit double linear fit up to, and over set radial value (r_fit, L30) then plot
#
coeffs = [0]*len(snaparr_tmp)
matcov = [0]*len(snaparr_tmp)
y_fitted = [[0]*len(sig[0,r_fit:r_limit] ) ]*len(snaparr_tmp)
#
for i in range(0, len(snaparr_tmp)):
#
	plt.figure(1)
	ax1 = plt.subplot(111)
#
	coeffs[i], matcov[i] = curve_fit(func, np.log10(r[snaparr_tmp[i],r_fit:r_limit]), \
	   np.log10(sig[snaparr_tmp[i],r_fit:r_limit]), [1, 1] )
	y_fitted[i] = func(np.log10(r[snaparr_tmp[i],r_fit:r_limit]), coeffs[i][0], coeffs[i][1])
#
	if (spline == "TRUE"):
		rnew = np.arange( r_start, r_limit, r_limit/(r_limit/spline_smooth) )
#
		tck = interpolate.splrep(r[snaparr_tmp[i],r_start:r_limit], \
		   sig[snaparr_tmp[i],r_start:r_limit])
		signew = interpolate.splev(rnew, tck, der = 0)
#
		line1 = plt.plot(np.log10(rnew), np.log10(signew), \
		   color = col_arr[i])
#
	elif (spline == "FALSE"):
		line1 = plt.plot(np.log10(r[snaparr_tmp[i],r_start:r_limit]), \
		   np.log10(sig[snaparr_tmp[i],r_start:r_limit]), \
		   color = col_arr[i])
#
	line2 = plt.plot(np.log10(r[snaparr_tmp[i],r_fit:r_limit]), y_fitted[i], \
	   color = col_arr[i], linewidth = 2, linestyle = 'dashed', \
	   label =str(round(coeffs[i][0], 4)) )
#
	plt.xlabel("log Radius", fontsize = 8) 
	plt.ylabel("log Surface density", fontsize = 8, labelpad=0.5)
	plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1])
	plt.yticks(fontsize = 8)
	plt.xticks(fontsize = 8)
	plt.legend(loc='upper right', title = (r'${\alpha}$')+" values of fit", fontsize=6)
#
	plt.savefig(str(plotdir)+'SD_r_'+str(r_limit)+'AU__DE05_'+str(snaparr[i]+69)+'.png')
	plt.clf()	
#
# Print power law fits to surface density profiles
#
if (powerfit_check == "TRUE"):
	print "\n Surface density power law fits \n"
	for i in range(0, len(snaparr_tmp)):
		print "Power index for the "+str(snaparr_tmp[i])+"th snapshot of " \
		   + str(ea_run)+" is: ",  \
		   str( round(coeffs[i][0] , 4) ) + " \n"
#
# Plot temperature vs. radius for singular accretion event
# First, fit double linear fit up to, and over set radial value (r_fit, L30) then plot
#
coeffs = [0]*len(snaparr_tmp)
matcov = [0]*len(snaparr_tmp)
y_fitted1 = [[0]*len(T[0,r_fit:r_limit] ) ]*len(snaparr_tmp)
#
for i in range(0, len(snaparr_tmp)):
#
	plt.figure(1)
	ax1 = plt.subplot(111)
#
	coeffs[i], matcov[i] = curve_fit(func, np.log10(r[snaparr_tmp[i],r_fit:r_limit]), \
	   np.log10(T[snaparr_tmp[i],r_fit:r_limit]), [1, 1])
	y_fitted1[i] = func(np.log10(r[snaparr_tmp[i],r_fit:r_limit]), coeffs[i][0], coeffs[i][1])
#
	if (spline == "TRUE"):
		rnew = np.arange(r_start,r_limit,r_limit/(r_limit/spline_smooth))
#
		tck = interpolate.splrep(r[snaparr_tmp[i],r_start:r_limit], \
		   T[snaparr_tmp[i],r_start:r_limit])
		Tnew = interpolate.splev(rnew, tck, der = 0)
#
		line1 = plt.plot(np.log10(rnew), np.log10(Tnew), \
		   color = col_arr[i])
#
	elif (spline == "FALSE"):
		line1 = plt.plot(np.log10(r[snaparr_tmp[i],r_start:r_limit]), \
		   np.log10(T[snaparr_tmp[i],r_start:r_limit]), color = col_arr[i])
#
	line2 = plt.plot(np.log10(r[snaparr_tmp[i],r_fit:r_limit]), y_fitted1[i], \
	   color = col_arr[i], linewidth = 2, linestyle = 'dashed', \
	   label = str(round(coeffs[i][0], 4)) )
#
	plt.xlabel("Disc radius (AU)", fontsize = 8) 
	plt.ylabel("log Temperature", fontsize = 8, labelpad=0.5)
	plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1])
	plt.yticks(fontsize = 8)
	plt.xticks(fontsize = 8)
	plt.legend(loc='upper right', title = (r'${\alpha}$')+" values of fit", fontsize=6)
#
	plt.savefig(str(plotdir)+'T_r_'+str(r_limit)+'AU__DE05_'+str(snaparr[i]+69)+'.png')
	plt.clf()	
#
# Print power law fits to Temperature profiles
#
if (powerfit_check == "TRUE"):
	print "\n Temperature power law fits \n"
	for i in range(0, len(snaparr_tmp)):
		print "Power index for the "+str(snaparr_tmp[i])+"th snapshot of "+ \
		   str(ea_run)+" is: ",  \
		   str( round(coeffs[i][0], 4) ) + " \n"
#
#
# Plot Q vs. radius
# First, fit double linear fit up to, and over set radial value (r_fit, L30) then plot
#
coeffs = [0]*len(snaparr_tmp)
matcov = [0]*len(snaparr_tmp)
y_fitted1 = [[0]*len(Q[0,r_fit:r_limit] ) ]*len(snaparr_tmp)
#
for i in range(0, len(snaparr_tmp)):
#
	plt.figure(1)
	ax1 = plt.subplot(111)
#
	coeffs[i], matcov[i] = curve_fit(func, r[snaparr_tmp[i],r_fit:r_limit], \
	   np.log10(Q[snaparr_tmp[i],r_fit:r_limit]), [1, 1])
	y_fitted1[i] = func(r[snaparr_tmp[i],r_fit:r_limit], coeffs[i][0], coeffs[i][1])
#
	if (spline == "TRUE"):
		rnew = np.arange(r_start,r_limit,r_limit/(r_limit/spline_smooth))
#
		tck = interpolate.splrep(r[snaparr_tmp[i],r_start:r_limit], \
		   Q[snaparr_tmp[i],r_start:r_limit])
		Qnew = interpolate.splev(rnew, tck, der = 0)
#
		line1 = plt.plot(np.log10(rnew), np.log10(Qnew), \
		   color = col_arr[i])
#
	elif (spline == "FALSE"):
		line1 = plt.plot(np.log10(r[snaparr_tmp[i],r_start:r_limit]), \
		   np.log10(Q[snaparr_tmp[i],r_start:r_limit]), color = col_arr[i])
#
	line2 = plt.plot(np.log10(r[snaparr_tmp[i],r_fit:r_limit]), y_fitted1[i], \
	   color = col_arr[i], linewidth = 2, linestyle = 'dashed', \
	   label = str(round(coeffs[i][0], 4)) )
#
	plt.xlabel("Disc radius (AU)", fontsize = 8) 
	plt.ylabel("log Q", fontsize = 8, labelpad=0.5)
	plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1])
	plt.yticks(fontsize = 8)
	plt.xticks(fontsize = 8)
	plt.legend(loc='upper right', title = (r'${\alpha}$')+" values of fit", fontsize=6)
#
	plt.savefig(str(plotdir)+'Q_r_'+str(r_limit)+'AU__DE05_'+str(snaparr[i]+69)+'.png')
	plt.clf()	
#
# Print power law fits to Temperature profiles
#
if (powerfit_check == "TRUE"):
	print "\n Q power law fits \n"
	for i in range(0, len(snaparr_tmp)):
		print "Power index for the "+str(snaparr_tmp[i])+"th snapshot of "+ \
		   str(ea_run)+" is: ",  \
		   str( round(coeffs[i][0], 4) ) + " \n"
#
#os.system('for f in '+plotdir+'*.pdf; do pdftops -eps $f; done')
#os.system('rm -r '+plotdir+'*.pdf')
