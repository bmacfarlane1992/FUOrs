#
# pdisc_r.py
#
# Python program to read in  pdisc (radially varied, temporally static) data structures. 
# Uses information of file name to generate plots analogous to Stamatellos, 
# Whitworth & Hubber (2012) work.
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
	# !!! DO TIME CHECK PRIOR TO ENTERING VALUES OF [tag]_s_[time] INDICES !!!
#
time_check = "FALSE"			# Choose whether ("TRUE") or not ("FALSE") to output times of accretion
powerfit_check = "FALSE"		# Choose whether ("TRUE") or not ("FALSE") to fitted power indices to Sig and T profiles
#
col_arr = ["b", "g", "r", "c", "m", "k"]
#
def func(x, a, b):			# Linear function used to fit power indices to log-log Sig/T vs. radius distributions
  return a*x + b
r_fit = 50				# Radial point after which power indices fit to data as described above
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
def read(arch_dir, plotdir, ea_run, snaparr, EA_timeref, EA_lenref, pmass, \
   pradius, hasharr_app, n_accr, r_limit, spline, r_d_kep, timearr, v_K):
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
		fcount = 0
		for a in range(0, len(snaparr)):
			if (snaparr[a] < (1000-70)):
				file_list.append(arch_dir+'pdisc/DE05.du.00'+ \
				   str(snaparr[a]+69)+'.pdisc.1')
			elif (snaparr[a] > (1000-70)):
				file_list.append(arch_dir+'pdisc/DE05.du.0'+ \
				   str(snaparr[a]+69)+'.pdisc.1')
			snaparr_tmp[a] = fcount
			fcount = fcount + 1
#
	elif (snaparr.ndim == 2):
		snaparr_tmp = np.array([[0]*len(snaparr[0])]*len(snaparr)) 
		fcount = 0
		for a in range(0, len(snaparr)):
			for b in range(0, len(snaparr[0])):
				if (snaparr[a][b] < (1000-70)):
					file_list.append(arch_dir+'pdisc/DE05.du.00'+ \
					   str(snaparr[a][b]+69)+'.pdisc.1')
				elif (snaparr[a][b] > (1000-70)):
					file_list.append(arch_dir+'pdisc/DE05.du.0'+ \
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
	time = np.array(time) ; r = np.array(r) ; Q = np.array(Q) ; T = np.array(T) ; sig = np.array(sig)
	rInM = np.array(rInM) ; M_s = np.array(M_s) ; M_d = np.array(M_d) ; v_r = np.array(v_r)
	v_the = np.array(v_the) ; v_z = np.array(v_z) ; vkep = np.array(vkep) ; v_mod = np.array(v_mod)
#
#
	# Plot the disc parameters before, during and after the accretion for defined radial value
#
	print("pdisc plots now being generated")
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# - - - NON-EA RUNS - - - # 
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
			line1 = plt.plot(r[snaparr_tmp[i],0:r_limit], Q[snaparr_tmp[i],0:r_limit], \
			   color = col_arr[i])
		plt.xlabel("Disc radius (AU)", fontsize = 8) 
		plt.ylabel("Toomre Q parameter", fontsize = 8, labelpad=0.5)
		plt.ylim(ax1.get_ylim()[0], 20)
		plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
#
		ax2 = plt.subplot(222)
		for i in range(0, len(snaparr_tmp)):
			line3 = plt.plot(r[snaparr_tmp[i],0:r_limit], rInM[snaparr_tmp[i],0:r_limit], \
			   color = col_arr[i])
		plt.xlabel("Disc radius (AU)", fontsize = 8) 
		plt.ylabel("Disc mass "+(r'(M$_{\odot}$)'), fontsize = 8, labelpad=0.5)
		plt.ylim(ax2.get_ylim()[0], ax2.get_ylim()[1])
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
#
		ax3 = plt.subplot(223)
		for i in range(0, len(snaparr_tmp)):
			line4 = plt.plot(r[snaparr_tmp[i],0:r_limit], v_r[snaparr_tmp[i],0:r_limit], \
			   color = col_arr[i])
		plt.xlabel("Disc radius (AU)", fontsize = 8) 
		plt.ylabel((r'v$_{r}$')+' (km'+(r's$^{-1}$')+')', fontsize = 8, labelpad=0.5)
		plt.ylim(ax3.get_ylim()[0], ax3.get_ylim()[1])
		plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
#
		plt.savefig(str(plotdir)+'Matrix_'+str(r_limit)+'AU'+'.pdf') ; plt.clf()
#
# Plot surface density vs. radius for singular accretion event
# First, fit double linear fit up to, and over set radial value (r_fit, L30) then plot
#
		coeffs1 = [0]*len(snaparr_tmp) ; coeffs2 = [0]*len(snaparr_tmp)
		matcov1 = [0]*len(snaparr_tmp) ; matcov2 = [0]*len(snaparr_tmp)
		y_fitted1 = [[0]*len(sig[0,0:r_limit] ) ]*len(snaparr_tmp)
		y_fitted2 = [[0]*len(sig[0,0:r_limit] ) ]*len(snaparr_tmp)
#
		plt.figure(1)
		ax1 = plt.subplot(111)
		for i in range(0, len(snaparr_tmp)):
#
			coeffs1[i], matcov1[i] = curve_fit(func, r[snaparr_tmp[i],5:r_fit], \
			   np.log10(sig[snaparr_tmp[i],5:r_fit]), [1, 1])
			y_fitted1[i] = func(r[snaparr_tmp[i],5:r_fit], coeffs1[i][0], coeffs1[i][1])
#
			coeffs2[i], matcov2[i] = curve_fit(func, r[snaparr_tmp[i],r_fit:r_limit], \
			   np.log10(sig[snaparr_tmp[i],r_fit:r_limit]), [1, 1])
			y_fitted2[i] = func(r[snaparr_tmp[i],r_fit:r_limit], coeffs2[i][0], \
			   coeffs2[i][1])
#
			if (spline == "TRUE"):
				rnew = np.arange(5,r_limit,r_limit/(r_limit/3))
#
				tck = interpolate.splrep(r[snaparr_tmp[i],5:r_limit], \
				   sig[snaparr_tmp[i],5:r_limit])
				signew = interpolate.splev(rnew, tck, der = 0)
#
				line1 = plt.plot(rnew, np.log10(signew), \
				   color = col_arr[i])
#
			elif (spline == "FALSE"):
				line1 = plt.plot(r[snaparr_tmp[i],0:r_limit], \
				   np.log10(sig[snaparr_tmp[i],0:r_limit]), \
			   	color = col_arr[i])
#
			line2 = plt.plot(r[snaparr_tmp[i],5:r_fit], y_fitted1[i], \
			   color = col_arr[i], linewidth = 2, linestyle = 'dashed', \
			   label = "-"+str(round(np.log10(abs(coeffs1[i][0])), 4)) )
			line3 = plt.plot(r[snaparr_tmp[i],r_fit:r_limit], y_fitted2[i], \
			   color = col_arr[i], linewidth = 2, linestyle = 'dashed', \
			   label = "-"+str(round(np.log10(abs(coeffs2[i][0])), 4)) )
#
		plt.xlabel("Disc radius (AU)", fontsize = 8) 
		plt.ylabel("log Surface density", fontsize = 8, labelpad=0.5)
		plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1])
		plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
		plt.legend(loc='upper right', title = (r'${\alpha}$')+" values of fit", fontsize=6)
#
		plt.savefig(str(plotdir)+'SD_r_'+str(r_limit)+'AU.pdf') ; plt.clf()	
#
# Print power law fits to surface density profiles
#
		if (powerfit_check == "TRUE"):
			print "\n Surface density power law fits \n"
			for i in range(0, len(snaparr_tmp)):
				print "Power index for the "+snaparr_tmp[i]+"th snapshot of " \
				   +str(ea_run)+" is: ",  \
				   str(round(np.log10(abs(coeffs1[i][0] )), 4)), \
				   str(round(np.log10(abs(coeffs2[i][0])), 4))+" \n"
#
# Plot temperature vs. radius for singular accretion event
# First, fit double linear fit up to, and over set radial value (r_fit, L30) then plot
#
		coeffs1 = [0]*len(snaparr_tmp) ; coeffs2 = [0]*len(snaparr_tmp)
		matcov1 = [0]*len(snaparr_tmp) ; matcov2 = [0]*len(snaparr_tmp)
		y_fitted1 = [[0]*len(T[0,0:r_limit] ) ]*len(snaparr_tmp)
		y_fitted2 = [[0]*len(T[0,0:r_limit] ) ]*len(snaparr_tmp)
#
		plt.figure(1)
		ax1 = plt.subplot(111)
		for i in range(0, len(snaparr_tmp)):
#
			coeffs1[i], matcov1[i] = curve_fit(func, r[snaparr_tmp[i],5:r_fit], \
			   np.log10(T[snaparr_tmp[i],5:r_fit]), [1, 1])
			y_fitted1[i] = func(r[snaparr_tmp[i],5:r_fit], coeffs1[i][0], coeffs1[i][1])
#
			coeffs2[i], matcov2[i] = curve_fit(func, r[snaparr_tmp[i],r_fit:r_limit], \
			   np.log10(T[snaparr_tmp[i],r_fit:r_limit]), [1, 1])
			y_fitted2[i] = func(r[snaparr_tmp[i],r_fit:r_limit], coeffs2[i][0], \
			   coeffs2[i][1])
			if (spline == "TRUE"):
				rnew = np.arange(5,r_limit,r_limit/(r_limit/3))
#
				tck = interpolate.splrep(r[snaparr_tmp[i],5:r_limit], \
				   T[snaparr_tmp[i],5:r_limit])
				Tnew = interpolate.splev(rnew, tck, der = 0)
#
				line1 = plt.plot(rnew, np.log10(Tnew), \
				   color = col_arr[i])
#
			elif (spline == "FALSE"):
				line1 = plt.plot(r[snaparr_tmp[i],0:r_limit], \
				   np.log10(T[snaparr_tmp[i],0:r_limit]), color = col_arr[i])
			line2 = plt.plot(r[snaparr_tmp[i],5:r_fit], y_fitted1[i], \
			   color = col_arr[i], linewidth = 2, linestyle = 'dashed', \
			   label = "-"+str(round(np.log10(abs(coeffs1[i][0])), 4)) )
			line3 = plt.plot(r[snaparr_tmp[i],r_fit:r_limit], y_fitted2[i], \
			   color = col_arr[i], linewidth = 2, linestyle = 'dashed', \
			   label = "-"+str(round(np.log10(abs(coeffs2[i][0])), 4)) )
#
		plt.xlabel("Disc radius (AU)", fontsize = 8) 
		plt.ylabel("log Temperature", fontsize = 8, labelpad=0.5)
		plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1])
		plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
		plt.legend(loc='upper right', title = (r'${\alpha}$')+" values of fit", fontsize=6)
#
		plt.savefig(str(plotdir)+'T_r_'+str(r_limit)+'AU.pdf') ; plt.clf()	
#
# Print power law fits to Temperature profiles
#
		if (powerfit_check == "TRUE"):
			print "\n Temperature power law fits \n"
			for i in range(0, len(snaparr_tmp)):
				print "Power index for the "+snaparr_tmp[i]+"th snapshot of "+ \
				   str(ea_run)+" is: ",  \
				   str(round(np.log10(abs(coeffs1[i][0])), 4)), \
				   str(round(np.log10(abs(coeffs2[i][0])), 4))+" \n"
#
# Plot azimuthal velocity vs. radius for singular accretion event. Compare data of before/during/
# after event with analytic calculation of Keplerian velocity for snapshot
#
		line1 = [ []*len(v_the[0,0:r_limit]) ]*len(snaparr_tmp)
		line2 = [ []*len(v_the[0,0:r_limit]) ]*len(snaparr_tmp)
		plt.figure(1)
		ax1 = plt.subplot(111)
		plt.title('Azimuthal velocity vs. radius')
		for i in range(0, len(snaparr_tmp)):
			line1[i] = plt.plot(r[snaparr_tmp[i],0:r_limit], \
			   abs(v_the[snaparr_tmp[i],0:r_limit]), color = col_arr[i])
			line2[i] = plt.plot(r[snaparr_tmp[i],0:r_limit], \
			   vkep[snaparr_tmp[i],0:r_limit], \
			   label = 'Keplerian velocity ('+str(snaparr_tmp[i])+' analytical)', \
			   linewidth = 2, linestyle = 'dashed', color=col_arr[i])
			plt.axvline(x=r_d_kep[snaparr_tmp[i]], linewidth = 2, \
			   linestyle = 'dotted', color=col_arr[i])
#
	# If a planet is present, overplot onto x-axis with hashed line
#
		if (len(pmass) > 0):
			for a in range(0, len(pmass)):
				for idash in range(0, len(pmass[0])):
					if (pmass[a][idash] != 0.):
						ax1.plot(pradius[a][idash], \
						   vkep[snaparr_tmp[idash],int(pradius[a][idash])], \
						    col_arr[idash]+'o', markersize=10.0)
		plt.xlabel("Disc radius (AU)", fontsize = 8)
		plt.ylabel(''+(r'v$_{\phi}$')+' (km'+(r's$^{-1}$')+')', fontsize = 8, labelpad=0.5)
		plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1]/2.) ; plt.xlim(0, r_limit)
		plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
		plt.legend(loc='upper right', fontsize=6)
#
		plt.savefig(str(plotdir)+'Azi_r_'+str(r_limit)+'AU.pdf') ; plt.clf()
#
	# Plot absolute velocity vs. keplerian expectation at radius R
#
# First, make sure arrays defined where velocity condition not met to create continuous fill plots
#
#
		n_ax = 311
		for a in range(0, len(snaparr_tmp)):
			ax1 = plt.subplot(n_ax)
			vcondarr1 = [] ; vcondarr2 = [] ; vcondarr3 = []
			jstart = 0
			for i in range(0, 401):
				for j in range(jstart, 401):
					if (abs(v_the[snaparr_tmp[a],j]) < \
					   (float(v_K)/100.)*abs(vkep[snaparr_tmp[a],j])):
						vcondarr1.append(j)
						for k in range(j, 401):
							if (abs(v_the[snaparr_tmp[a],k]) \
							   > (float(v_K)/100.)*abs(vkep[snaparr_tmp[a],k])):
								vcondarr1.append(k) ; jstart = k ; break
						break
			vcondarr1_app=[] ; vcondarr2_app=[] ; vcondarr3_app=[]
			for i in vcondarr1:
	 			if i not in vcondarr1_app:
	       				vcondarr1_app.append(i)
			if (len(vcondarr1_app) % 2 == 1):
				vcondarr1_app.append(401)
			n_cond1 = len(vcondarr1_app)/2	
			vcondarr1_app = np.reshape(vcondarr1_app, (n_cond1, 2))	
#	
# Now begin plotting
#
			line1 = plt.plot(r[snaparr_tmp[a],0:r_limit], \
			   abs(v_the[snaparr_tmp[a],0:r_limit]), label="Azimuthal (data)", \
			   color = col_arr[0])
			line11 = plt.plot(r[snaparr_tmp[a],0:r_limit], \
			   vkep[snaparr_tmp[a],0:r_limit], \
			   label = 'Keplerian velocity (analytical)', \
			   linewidth = 2, linestyle = 'dashed', color = col_arr[1])
			line111 = plt.plot(r[snaparr_tmp[a],0:r_limit], \
			   v_mod[snaparr_tmp[a],0:r_limit], label="Mag. velocity (data)", \
			   color = col_arr[2])
			line1111 = plt.plot(r[snaparr_tmp[a],0:r_limit], \
			   abs(v_r[snaparr_tmp[a],0:r_limit]), label="Radial velocity (data)", \
			   color = col_arr[3])
			line11111 = plt.plot(r[snaparr_tmp[a],0:r_limit], \
			   abs(v_z[snaparr_tmp[a],0:r_limit]), label="Z velocity (data)", \
			   color = col_arr[4])
			for i in range(0,n_cond1):
				plt.fill_between(r[snaparr_tmp[a], \
				   vcondarr1_app[i][0]:vcondarr1_app[i][1]], \
				   ax1.get_ylim()[0], ax1.get_ylim()[1], \
				   color=col_arr[5], alpha = 0.5)
			plt.xlabel("Disc radius (AU)") ; plt.ylabel("|v|"+' (km'+(r's$^{-1}$')+')')
			plt.xlim(0,r_limit) ; plt.ylim(0., max(vkep[snaparr_tmp[a],1:r_limit]))
			plt.title(str(round(timearr[a], 2))+") kyr snapshot")
			plt.legend(loc='upper right', fontsize=12)
			n_ax = n_ax+1
#
		plt.savefig(str(plotdir)+'vthe_restr.pdf') ; plt.clf()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# - - - EA RUNS [ VARYING OUTBURST DURATION ]- - - # 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# For EA runs, plot parameters comparing SHORT, MEDIUM and LONG snapshots
#
	if (snaparr_tmp.ndim == 2):
#
	# Check the times of snapshots entered for before, during and after accretion event chosen
#
		if (time_check == "TRUE"):
			print "\n"
			for i in range(0,n_accr):
				print "Accretion event "+str(i+1)+" begins at snapshot "+ \
				   str(hasharr_app[i][0])+" for "+ \
				   str((hasharr_app[i][1]-hasharr_app[i][0])*10)+" years" \
				   " until snapshot "+str(hasharr_app[i][1])
			print "\n The final snapshot of ea run "+str(ea_run)+" is "+str(len(time))+"\n"
			for i in range(0, len(snaparr_tmp)):
				for j in range(0, len(snaparr_tmp[0])):
					print "The time of snapshot ("+EA_lenref[i]+", "+ \
					   EA_timeref[j]+") entered is: ", round(timearr[i][j], 3)
			print "\n"
			exit()
#
	# Invoke a switch statement to loop over pdisc analysis of ea [2, 3, 4, 5, 6] runs
	# with focus on (a) single accretion event and (b) focus on snapshot relative to accretion
	# state for different events
#		
		for lendur in range(0, 2):
			title_point = [""]*3
			leg_point = [""]*3
			if (lendur == 0):
				title_point = EA_lenref
				leg_point = EA_timeref
			elif (lendur == 1):
				title_point = EA_timeref
				title_point = EA_lenref
				snaparr_tmp = np.rot90(snaparr_tmp, 2)
#
		for i in range(0,len(snaparr_tmp)):
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
			plt.xlabel("Disc radius (AU)", fontsize = 8) 
			plt.ylabel("Toomre Q parameter", fontsize = 8, labelpad=0.5)
			plt.ylim(ax1.get_ylim()[0], 20)
			plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
#
			ax2 = plt.subplot(222)
			for j in range(0, len(snaparr_tmp[0])):
				line2 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
				   np.log10(T[snaparr_tmp[i][j],0:r_limit]), \
				   label = "t = "+str(round(timearr[i][j], 2))+" kyr "+ \
				   EA_timeref[j], color = col_arr[j])
			plt.xlabel("Disc radius (AU)", fontsize = 8) 
			plt.ylabel("log Temperature", fontsize = 8, labelpad=0.5)
			plt.ylim(ax2.get_ylim()[0], ax2.get_ylim()[1])
 			plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
			plt.legend(loc='upper right', fontsize=6)
#
			ax3 = plt.subplot(223)
			for j in range(0, len(snaparr_tmp[0])):
				line3 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
				   rInM[snaparr_tmp[i][j],0:r_limit], color = col_arr[j])
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.xlabel("Disc radius (AU)", fontsize = 8) 
			plt.ylabel("Disc mass "+(r'(M$_{\odot}$)'), fontsize = 8, labelpad=0.5)
			plt.ylim(ax3.get_ylim()[0], ax3.get_ylim()[1])
			plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
#
			ax4 = plt.subplot(224)
			for j in range(0, len(snaparr_tmp[0])):
				line4 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
				   v_r[snaparr_tmp[i][j],0:r_limit], color = col_arr[j])
			plt.xlabel("Disc radius (AU)", fontsize = 8) 
			plt.ylabel((r'v$_{r}$')+' (km'+(r's$^{-1}$')+')',fontsize = 8,labelpad=0.5)
			plt.ylim(ax4.get_ylim()[0], ax4.get_ylim()[1])
			plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
#
			plt.savefig(str(plotdir)+'Matrix_'+str(r_limit)+'AU_'+EA_lenref[i]+'.pdf')
			plt.clf()
#
# Plot surface density vs. radius for comparative accretion reference times. 
# First, fit double linear fit up to, and over set radial value (r_fit, L30) then plot
#
			coeffs1 = [0]*len(snaparr_tmp[0]) ; coeffs2 = [0]*len(snaparr_tmp[0])
			matcov1 = [0]*len(snaparr_tmp[0]) ; matcov2 = [0]*len(snaparr_tmp[0])
			y_fitted1 = [[0]*len(sig[snaparr_tmp[0][0],5:r_limit] ) ]*len(snaparr_tmp[0])
			y_fitted2 = [[0]*len(sig[snaparr_tmp[0][0],r_fit:r_limit] ) ]*len(snaparr_tmp[0])
#
			plt.figure(1)
#
			ax1 = plt.subplot(111)
#
			for j in range(0, len(snaparr_tmp[0])):
				coeffs1[j], matcov1[j] = curve_fit(func, r[snaparr_tmp[i][j],5:r_fit], \
				   np.log10(sig[snaparr_tmp[i][j],5:r_fit]), [1, 1])
				y_fitted1[j] = func(r[snaparr_tmp[i][j],5:r_fit], \
				   coeffs1[j][0], coeffs1[j][1])
#	
				coeffs2[j], matcov2[j] = curve_fit(func, r[snaparr_tmp[i][j],r_fit:r_limit], \
				   np.log10(sig[snaparr_tmp[i][j],r_fit:r_limit]), [1, 1])
				y_fitted2[j] = func(r[snaparr_tmp[i][j],r_fit:r_limit], \
				   coeffs2[j][0], coeffs2[j][1])
				#
				if (spline == "TRUE"):
					rnew = np.arange(5,r_limit,r_limit/(r_limit/3))
#
					tck = interpolate.splrep(r[snaparr_tmp[i][j],5:r_limit], \
					   sig[snaparr_tmp[i][j],5:r_limit])
					signew = interpolate.splev(rnew, tck, der = 0)
#
					line1 = plt.plot(rnew, np.log10(signew), \
					   color = col_arr[j])
#
				elif (spline == "FALSE"):
					line1 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
					   np.log10(sig[snaparr_tmp[i][j],0:r_limit]), \
				   	   color = col_arr[j])
#
				line2 = plt.plot(r[snaparr_tmp[i][j],5:r_fit], y_fitted1[j], 
				   color = col_arr[j], linewidth = 2, linestyle = 'dashed', \
				   label = "-"+str(round(np.log10(abs(coeffs1[j][0])), 4)) )
				line3 = plt.plot(r[snaparr_tmp[i][j],r_fit:r_limit], y_fitted2[j], \
				   color = col_arr[j], linewidth = 2, linestyle = 'dashed', \
				   label = "-"+str(round(np.log10(abs(coeffs2[j][0])), 4)) )
#
			plt.xlabel("Disc radius (AU)", fontsize = 8) 
			plt.ylabel("log Surface Density", fontsize = 8, labelpad=0.5)
			plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1])
			plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
			plt.legend(loc='upper right', title = (r'${\alpha}$')+" values of fit", \
			   fontsize=6)
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
					   +"th snapshot of "+str(ea_run)+" is: ",  \
					   str(round(np.log10(abs(coeffs1[j][0] )), 4)), \
					   str(round(np.log10(abs(coeffs2[j][0])), 4))+" \n"
#
# Plot temperature vs. radius for comparative accretion reference times. 
# First, fit double linear fit up to, and over set radial value (r_fit, L30) then plot
#
			coeffs1 = [0]*len(snaparr_tmp[0]) ; coeffs2 = [0]*len(snaparr_tmp[0])
			matcov1 = [0]*len(snaparr_tmp[0]) ; matcov2 = [0]*len(snaparr_tmp[0])
			y_fitted1 = [[0]*len(T[snaparr_tmp[0][0],0:r_limit] ) ]*len(snaparr_tmp[0])
			y_fitted2 = [[0]*len(T[snaparr_tmp[0][0],0:r_limit] ) ]*len(snaparr_tmp[0])
#
			plt.figure(1)
#
			ax1 = plt.subplot(111)
#
			for j in range(0, len(snaparr_tmp[0])):
				coeffs1[j], matcov1[j] = curve_fit(func, r[snaparr_tmp[i][j],5:r_fit], \
				   np.log10(T[snaparr_tmp[i][j],5:r_fit]), [1, 1])
				y_fitted1[j] = func(r[snaparr_tmp[i][j],5:r_fit], \
				   coeffs1[j][0], coeffs1[j][1])
#
				coeffs2[j], matcov2[j] = curve_fit(func, r[snaparr_tmp[i][j],r_fit:r_limit], \
				   np.log10(T[snaparr_tmp[i][j],r_fit:r_limit]), [1, 1])
				y_fitted2[j] = func(r[snaparr_tmp[i][j],r_fit:r_limit], \
				   coeffs2[j][0], coeffs2[j][1])
#
				if (spline == "TRUE"):
					rnew = np.arange(5,r_limit,r_limit/(r_limit/3))
#
					tck = interpolate.splrep(r[snaparr_tmp[i][j],5:r_limit], \
					   T[snaparr_tmp[i][j],5:r_limit])
					Tnew = interpolate.splev(rnew, tck, der = 0)
#
					line1 = plt.plot(rnew, np.log10(Tnew), \
					   color = col_arr[j])
#
				elif (spline == "FALSE"):
					line1 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
					   np.log10(T[snaparr_tmp[i][j],0:r_limit]), \
				   	   color = col_arr[j])
#
				line2 = plt.plot(r[snaparr_tmp[i][j],5:r_fit], y_fitted1[j], 
				   color = col_arr[j], linewidth = 2, linestyle = 'dashed', \
				   label = "-"+str(round(np.log10(abs(coeffs1[j][0])), 4)) )
				line3 = plt.plot(r[snaparr_tmp[i][j],r_fit:r_limit], y_fitted2[j], \
				   color = col_arr[j], linewidth = 2, linestyle = 'dashed', \
				   label = "-"+str(round(np.log10(abs(coeffs2[j][0])), 4)) )
#
			plt.xlabel("Disc radius (AU)", fontsize = 8) 
			plt.ylabel("log Temperature", fontsize = 8, labelpad=0.5)
			plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1])
			plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
			plt.legend(loc='upper right', title = (r'${\alpha}$')+" values of fit", \
			   fontsize=6)
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
					   + "th snapshot of "+str(ea_run)+" is: ",  \
					   str(round(np.log10(abs(coeffs1[j][0] )), 4)), \
					   str(round(np.log10(abs(coeffs2[j][0])), 4))+" \n"
#
# Compare velocity data with analytic calculation of Keplerian velocity
#
			plt.figure(1)
			ax1 = plt.subplot(111)
			plt.title('Azimuthal velocity vs. radius for '+str(title_point[i])+' run')
#
			for j in range(0, len(snaparr_tmp[0])):

				line1 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], \
				   abs(v_the[snaparr_tmp[i][j],0:r_limit]), color=col_arr[j])
				line2 = plt.plot(r[snaparr_tmp[i][j],0:r_limit], 
				   vkep[snaparr_tmp[i][j],0:r_limit], \
				   label = 'Keplerian velocity ('+str(leg_point[j])+' analytical)', \
				   linewidth = 2, linestyle = 'dashed', color=col_arr[j])
#
	# If a planet is present, overplot onto x-axis with hashed line
#
				if (len(pmass) > 0):
					for a in range(0, len(pmass)):
						if (pmass[a][i][j] != 0.):
							ax1.plot(pradius[a][i][j],
							   vkep[snaparr_tmp[i][j],int(pradius[a][i][j])], \
							   col_arr[j]+'o', markersize=10.0)
#
				if (r_d_kep != 0):
					plt.axvline(x=r_d_kep[snaparr_tmp[i][j]], linewidth = 2, \
					   linestyle = 'dotted', color = col_arr[j])
#
			plt.xlabel("Disc radius (AU)", fontsize = 8)
			plt.ylabel(''+(r'v$_{\phi}$')+' (km'+(r's$^{-1}$')+')', \
			   fontsize = 8, labelpad=0.5)
			plt.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1]/2.) ; plt.xlim(0, r_limit)
			plt.yticks(fontsize = 8) ; plt.xticks(fontsize = 8)
			plt.legend(loc='upper right', fontsize=6)
#
			plt.savefig(str(plotdir)+'Azi_r_'+str(r_limit)+'AU_'+str(title_point[i])+'.pdf')
			plt.clf()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# - - - PLOT ABSOLUTE VELOCITY COMPARED TO EXPECTED KEPLERIAN VELOCITY PROFILES - - - # 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
# First, make sure arrays defined where velocity condition not met to create continuous fill plots
#
		for a in range(0, len(snaparr_tmp)):
			fig = plt.figure(2, figsize=(20,15))
#
			for b in range(0, len(snaparr_tmp[0])):
				ax1 = plt.subplot(311+b)
				vcondarr1 = []   ;   vcondarr2 = []   ;   vcondarr3 = []
				jstart = 0
				for i in range(0, 401):
					for j in range(jstart, 401):
						if (abs(v_the[snaparr_tmp[a][b],j]) < \
						   (float(v_K)/100.)*abs(vkep[snaparr_tmp[a][b],j])):
							vcondarr1.append(j)
							for k in range(j, 401):
								if (abs(v_the[snaparr_tmp[a][b],k]) \
								   > (float(v_K)/100.)*abs(vkep[snaparr_tmp[a][b],k])):
									vcondarr1.append(k)
								jstart = k ; break
							break
				vcondarr1_app=[] ; vcondarr2_app=[] ; vcondarr3_app=[]
				for i in vcondarr1:
		 			if i not in vcondarr1_app:
		       				vcondarr1_app.append(i)
				if (len(vcondarr1_app) % 2 == 1):
					vcondarr1_app.append(401)
				n_cond1 = len(vcondarr1_app)/2	
				vcondarr1_app = np.reshape(vcondarr1_app, (n_cond1, 2))	
#	
	# Now begin plotting
#
				line1 = plt.plot(r[snaparr_tmp[a][b],0:r_limit], \
				   abs(v_the[snaparr_tmp[a][b],0:r_limit]), \
				   label="Azimuthal (data)", color = col_arr[0])
				line11 = plt.plot(r[snaparr_tmp[a][b],0:r_limit], \
				   vkep[snaparr_tmp[a][b],0:r_limit], \
				   label = 'Keplerian velocity ('+str(EA_timeref[b])+' analytical)', \
				   linewidth = 2, linestyle = 'dashed', color = col_arr[1])
				line111 = plt.plot(r[snaparr_tmp[a][b],0:r_limit], \
				   v_mod[snaparr_tmp[a][b],0:r_limit], \
				   label="Mag. velocity (data)", color = col_arr[2])
				line1111 = plt.plot(r[snaparr_tmp[a][b],0:r_limit], \
				   abs(v_r[snaparr_tmp[a][b],0:r_limit]), \
				   label="Radial velocity (data)", color = col_arr[3])
				line11111 = plt.plot(r[snaparr_tmp[a][b],0:r_limit], \
				   abs(v_z[snaparr_tmp[a][b],0:r_limit]), \
				   label="Z velocity (data)", color = col_arr[4])
				for i in range(0,n_cond1):
					plt.fill_between(r[snaparr_tmp[a][b], \
					   vcondarr1_app[i][0]:vcondarr1_app[i][1]], \
					   ax1.get_ylim()[0], ax1.get_ylim()[1], \
					   color=col_arr[5], alpha = 0.5)
				plt.xlabel("Disc radius (AU)") ; plt.ylabel("|v|"+' (km'+(r's$^{-1}$')+')')
				plt.xlim(0,r_limit) ; plt.ylim(0., max(vkep[snaparr_tmp[a][b],1:r_limit]))
				plt.title(str(EA_timeref[b]) \
				   +" accretion event ("+str(round(timearr[a][b], 2))+") kyr)")
				plt.legend(loc='upper right', fontsize=12)
#
			plt.savefig(str(plotdir)+'vthe_restr_'+EA_lenref[a]+'.pdf') ; plt.clf()
#
	return r, vkep
