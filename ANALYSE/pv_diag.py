#
# pv_diag.py
#
# Programme designed to read in simulation data and generate synthetic PV diagram 
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
	# Unit conversion parameters
#
pcAU = 206265.
G = 6.67e-11
kg_Msol = 1.998e30
#
	# Miscellaneous parameters
#
kep_fit = "TRUE"
#
fit_start = 1
fit_stop = 10
#
RL_fitcheck = "FALSE"
#
raw_fit = "FALSE"
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
def pv(arch_dir, plotdir, ea_run, snaparr, v_K, inclin, r, vkep, EA_lenref, EA_timeref):
#
	print "Position-Velocity diagram now being plotted"
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Set filenames to be read in for PV plotting/analysis, dependent on EA run being analysed #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# For EA [0, 1] runs
#
	if (snaparr.ndim == 1):
		file_n = len(snaparr)
		pv_file = [""]*file_n
		Lfit_file = [""]*file_n
		Rfit_file = [""]*file_n
#
		snaparr_tmp = [0]*file_n
		plot_title = [""]*file_n ; plot_filename = [""]*file_n
		fcount = 0
		for i in range(0, len(snaparr)):
#
			if (snaparr[i] < (1000-70)):
				pv_file[fcount] = arch_dir+'pv_diag/DE05.du.00'+ \
				   str(snaparr[i]+69)+'.du.hist_pv.1'
				Lfit_file[fcount] = arch_dir+'pv_diag/DE05.du.00'+ \
				   str(snaparr[i]+69)+'.du.fit_Lpv.1'
				Rfit_file[fcount] = arch_dir+'pv_diag/DE05.du.00'+ \
				   str(snaparr[i]+69)+'.du.fit_Rpv.1'
			elif (snaparr[i] > (1000-70)):
				pv_file[fcount] = arch_dir+'pv_diag/DE05.du.0'+ \
				   str(snaparr[i]+69)+'.du.hist_pv.1'
				Lfit_file[fcount] = arch_dir+'pv_diag/DE05.du.0'+ \
				   str(snaparr[i]+69)+'.du.fit_Lpv.1'
				Rfit_file[fcount] = arch_dir+'pv_diag/DE05.du.0'+ \
				   str(snaparr[i]+69)+'.du.fit_Rpv.1'
#
			plot_title[fcount] = 'P-V Diagram: '+str(ea_run)+'   ('+str(snaparr[i])+')'
			plot_filename[fcount] = plotdir+'TEST_PV_diag_'+str(snaparr[i])+'.png'
#
			snaparr_tmp[fcount] = snaparr[i]
			fcount = fcount + 1
#
	# For EA [2, 3, 4, 5, 6] runs
#
	elif (snaparr.ndim == 2):
		file_n = len(snaparr)*len(snaparr[0])
		pv_file = [""]*file_n
		Lfit_file = [""]*file_n
		Rfit_file = [""]*file_n
#
		snaparr_tmp = [0]*file_n
		plot_title = [""]*file_n ; plot_filename = [""]*file_n
		fcount = 0
		for i in range(0, len(snaparr)):
			for j in range(0, len(snaparr[0])):
#
				if (snaparr[i][j] < (1000-70)):
					pv_file[fcount] = arch_dir+'pv_diag/DE05.du.00'+ \
					   str(snaparr[i][j]+69)+'.du.hist_pv.1'
					Lfit_file[fcount] = arch_dir+'pv_diag/DE05.du.00'+ \
					   str(snaparr[i][j]+69)+'.du.fit_Lpv.1'
					Rfit_file[fcount] = arch_dir+'pv_diag/DE05.du.00'+ \
					   str(snaparr[i][j]+69)+'.du.fit_Rpv.1'
				elif (snaparr[i][j] > (1000-70)):
					pv_file[fcount] = arch_dir+'pv_diag/DE05.du.0'+ \
					   str(snaparr[i][j]+69)+'.du.hist_pv.1'
					Lfit_file[fcount] = arch_dir+'pv_diag/DE05.du.0'+ \
					   str(snaparr[i][j]+69)+'.du.fit_Lpv.1'
					Rfit_file[fcount] = arch_dir+'pv_diag/DE05.du.0'+ \
					   str(snaparr[i][j]+69)+'.du.fit_Rpv.1'
#
				plot_title[fcount] = 'P-V Diagram: ('+str(EA_lenref[i])+', '+ \
				   str(EA_timeref[j])+')'
				plot_filename[fcount] = plotdir+'PV_diag_'+str(EA_lenref[i])+'_' \
				   +str(EA_timeref[j])+'.png'
#
				snaparr_tmp[fcount] = snaparr[i][j]
				fcount = fcount + 1
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Define arrays to fill for comparison of PV fitted masses to simulation criteria masses #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	pv_mass = []
	pv_kepcrit = []
	pv_pivcrit = []
	pv_sigALMAcrit = []
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Loop over list of snapshots and read in data #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	for i in range(0, file_n):
#
		r_bin = [] ; vz_bin = [] ; count_bin = []
#
	# Read data files, for PV diagram data, and data from which Keplerian fits are made
#	
		f = open(pv_file[i], 'r')
		r_range = 2. * float(f.readline()) ; v_range = 2. * float(f.readline())
		delt_r = 2. * float(f.readline()) ; delt_v = 2. * float(f.readline())
		y_min = 2. * float(f.readline()) ; y_max = 2. * float(f.readline())
		for line in f:
			line = line.strip() ; columns = line.split()
			r_bin.append(float(columns[0]))
			vz_bin.append(float(columns[1]))
			count_bin.append(float(columns[2]))
		f.close()
#
		r_bin = np.array(r_bin) ; vz_bin = np.array(vz_bin)
		count_bin = np.array(count_bin)
#
	# Select whether or not Keplerian fits are made to PV data
#
		if (kep_fit == "TRUE"):
#
	# Calculate the noise from which threshold of Keplerian profile is determined from
	# (see Seifried et al., 2016 for details)
#
			noise_count = 0
			for j in range(0, len(count_bin)):
				if (count_bin[j] != 0):
					noise_count = noise_count + 1
			sig_noise = math.sqrt( math.fsum(count_bin) / noise_count)
#	
			sig_noise = 1 * sig_noise
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# LHS PV diagram analysis #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
			Lfit_arr = [[],[],[]]
#
			f = open(Lfit_file[i], 'r')
			for line in f:
				line = line.strip() ; columns = line.split()
				Lfit_arr[0].append(float(columns[0]))
				Lfit_arr[1].append(float(columns[1]))
				Lfit_arr[2].append(float(columns[2]))
			f.close()
#	
			r_iter = int(math.ceil(2. * r_range / delt_r))
			v_iter = int(math.ceil(2. * v_range / delt_v))
#	
	# Define radial points of analysis, to begin evaluating y-positions of fits
#
			r_iter = r_iter/2 ; r_arr = [0.]*r_iter	
			rcount = 0
			for a in range(0, r_iter):
				rcount = rcount + 1
#
	# Determine the velocity value of Keplerian edge for fitting procedure
	# Ensure radial bins not at R = 0, as fitting routine crashes
#
				r_arr[a] = ((rcount-1)*delt_r) - r_range
				Lfit_y = []
				Lfit_r = []
#
				for a in range(0, r_iter):
					Lfit_ydum = []
					Lfit_ndum = []
					for b in range(0, len(Lfit_arr[0])):
						if (Lfit_arr[0][b] == r_arr[a]):
							Lfit_ydum.append(Lfit_arr[1][b])
							Lfit_ndum.append(Lfit_arr[2][b])
					for b in range(0, len(Lfit_ydum)):
						if (Lfit_ndum[b] > sig_noise):
							Lfit_y.append(Lfit_ydum[b])
							if (r_arr[a] == 0):
								Lfit_r.append(1e-1)
							else:
								Lfit_r.append(r_arr[a])
							break
#
	# Change quadrant of LHS fit arrays for correct power index fitting to curve
#
			Lfit_r = [-k for k in Lfit_r] ; Lfit_y = [-k for k in Lfit_y]
			Lfit_r = Lfit_r[::-1] ; Lfit_y = Lfit_y[::-1]
#
	# Define function of Keplerian fit
#	
			def func(x, a, b):
			  return a* np.power(x, -0.5) + b
#
	# Fitting routine to PV-diagram, and determination of system mass using Keplerian edge fit
#
			coeffs1, matcov1 = curve_fit(func, Lfit_r[fit_start:fit_stop], \
			   Lfit_y[fit_start:fit_stop], [1., 1.])
#
			Lfit_r1 = np.array(Lfit_r)*1.5e11 ; Lfit_y1 = np.array(Lfit_y)*1000.
			coeffs11, matcov11 = curve_fit(func, Lfit_r1[fit_start:fit_stop], \
			   Lfit_y1[fit_start:fit_stop], [1e10, 1e3])
			y_fitted1 = func(Lfit_r[fit_start:fit_stop], coeffs1[0], coeffs1[1])
			M_sys1 = (coeffs11[0]**2.0)/(G*kg_Msol)
#
	# Re-convert arrays to original quadrant, and change
#
			Lfit_r = [-k for k in Lfit_r] ; Lfit_y = [-k for k in Lfit_y]
			Lfit_r = Lfit_r[::-1] ; Lfit_y = Lfit_y[::-1]
#
			y_fitted1 = [-k for k in y_fitted1]
			y_fitted1 = y_fitted1[::-1]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# RHS PV diagram analysis #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
			Rfit_arr = [[],[],[]]
#
			f = open(Rfit_file[i], 'r')
			for line in f:
				line = line.strip() ; columns = line.split()
				Rfit_arr[0].append(float(columns[0]))
				Rfit_arr[1].append(float(columns[1]))
				Rfit_arr[2].append(float(columns[2]))
			f.close()
#
			r_iter = int(math.ceil(2. * r_range / delt_r))
			v_iter = int(math.ceil(2. * v_range / delt_v))
#
			r_iter = r_iter/2 ; r_arr = [0.]*r_iter	
			rcount = 0
			for a in range(0, r_iter):
				rcount = rcount + 1
#
				r_arr[a] = ((rcount-1)*delt_r)
				Rfit_y = []
				Rfit_r = []
#
				for a in range(0, r_iter):
					Rfit_ydum = []
					Rfit_ndum = []
					for b in range(0, len(Rfit_arr[0])):
						if (Rfit_arr[0][b] == r_arr[a]):
							Rfit_ydum.append(Rfit_arr[1][b])
							Rfit_ndum.append(Rfit_arr[2][b])
					ind = len(Rfit_ydum)-1
					for b in range(0, len(Rfit_ydum)):
						if (Rfit_ndum[ind] > sig_noise):
							Rfit_y.append(Rfit_ydum[ind])
							if (r_arr[a] == 0):
								Rfit_r.append(1e-1)
							else:
								Rfit_r.append(r_arr[a])
							break
						ind = ind - 1
#
			coeffs2, matcov2 = curve_fit(func, Rfit_r[fit_start:fit_stop], \
			   Rfit_y[fit_start:fit_stop], [1., 1.])
#
			Rfit_r2 = np.array(Rfit_r)*1.5e11 ; Rfit_y2 = np.array(Rfit_y)*1000.
			coeffs22, matcov22 = curve_fit(func, Rfit_r2[fit_start:fit_stop], \
			   Rfit_y2[fit_start:fit_stop], [1e10, 1e3])
			y_fitted2 = func(Rfit_r[fit_start:fit_stop], coeffs2[0], coeffs2[1])
			M_sys2 = (coeffs22[0]**2.0)/(G*kg_Msol)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# Outputs of PV diagram mass analysis #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
			pv_mass.append( (M_sys1 + M_sys2) / 2. )
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# Outputs of PV diagram mass analysis #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
			if (RL_fitcheck == "TRUE"):
				print "Mass of the system determined from LHS PV fit is: ", \
				   round(M_sys1,3), " Solar masses"
				print "Mass of the system determined from RHS PV fit is: ", \
				   round(M_sys2,3), " Solar masses"
				print "Averaged mass from PV diagram fits is: ", \
				   round(pv_mass[i],3), "Solar masses"
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# Plotting of PV diagram #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# Set up grid of interpolation points, interpolate on arbitrary grid, and plot
#
		xi, yi = np.linspace(r_bin.min(), r_bin.max(), 500), \
		   np.linspace(vz_bin.min(), vz_bin.max(), 500)
		xi, yi = np.meshgrid(xi, yi)
#
		zi = scipy.interpolate.griddata((r_bin, vz_bin), count_bin, (xi, yi), \
		   method='nearest')
#
		plt.figure(1)
		ax1 = plt.subplot(111)
		plt.imshow(zi, vmin=count_bin.min(), vmax=count_bin.max(), origin='lower', \
		   extent=[r_bin.min(), r_bin.max(), vz_bin.min(), vz_bin.max()], aspect = 'auto')
#
	# Overplot Keplerian fits to data, plotting points fitted to if chosen
#
		if ((raw_fit == "TRUE") and (kep_fit == "TRUE")):
			Lraw = plt.plot(Lfit_r, Lfit_y, linewidth = 2, \
			   linestyle = 'dashed', color = 'g')		
			Rraw = plt.plot(Rfit_r, Rfit_y, linewidth = 2, \
			   linestyle = 'dashed', color = 'g')		
#
		if (kep_fit == "TRUE"):
#
	# PV Keplerian distributions
#
			Lfit_start = len(Lfit_r)-fit_stop
			Lfit_stop = len(Lfit_r)-fit_start
			Lfit_start = len(Lfit_r)-fit_stop + (len(Lfit_r[Lfit_start:Lfit_stop]) - \
			   len(y_fitted1))
#
			Lfit = plt.plot(Lfit_r[Lfit_start:Lfit_stop], y_fitted1, linewidth = 2, \
			   linestyle = 'dashed', color = 'r')
			Rfit = plt.plot(Rfit_r[fit_start:fit_stop], y_fitted2, linewidth = 2, \
			   linestyle = 'dashed', color = 'r')		
#
	# Simulation Keplerian distributions
#
			plt.plot(r[i][3:75], vkep[i][3:75], linewidth = 2, \
			   linestyle = 'dashed', color = 'k')
			r[i] = [-k for k in r[i]] ; vkep[i] = [-k for k in vkep[i]]
			plt.plot(r[i][3:75], vkep[i][3:75], linewidth = 2, \
			   linestyle = 'dashed', color = 'k')
#
		plt.colorbar()
		plt.xlabel('Radius (AU)')
		plt.ylabel('Line of sight velocity'+' (km'+(r's$^{-1}$')+')')
		plt.xlim(-100,100)
		plt.ylim(-10, 10)
		plt.title(plot_title[i])
		plt.savefig(plot_filename[i])
		plt.clf()
#	
	return pv_mass
#
#
