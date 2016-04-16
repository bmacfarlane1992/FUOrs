#
# main.py
#
# Python program to execute read and analysis of simulation files.
# See script headers for more details on code I/O.
#
# Author: Benjamin MacFarlane
# Date: 16/04/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
arch_dir = "/home/ben/Documents/WORK_PLANETS/PROJECTS/FUORS/"	# Location of ea/ directory with data set
#
v_K = ["90"]			# Keplerian velocity percentage restriction on disc mass/radius
inclin = ["0","60","90"]		# Inclination of disc being analysed
#
ea_run = [0,1,3]		# Select EA runs to process
#
r_limit = 150			# Limit of radial plots in pdisc analyses
#
r_inspec = 50 			# Radius at which disc parameters inspected in rdisc
#
pv = "TRUE"			# Choose whether ("TRUE") or not ("FALSE") to generate PV diagram
mcomp_tmp = "TRUE"		# Choose whether ("TRUE") or not ("FALSE") to generate mass comparison
#				# of simulation vs. PV diagram analysis system masses
exp = "TRUE"			# Choose whether ("TRUE") or not ("FALSE") to compare disc radial SD and T profile exponents

#
EA_timeref = ["BEFORE","DURING_A","DURING_B","AFTER"]		# Define names of EA snapshots for EA length reference	
EA_lenref = ["SHORT","MEDIUM","LONG"]				# and time reference to EA outburst event
#
d = {}
snap0 = 'snap0' ; snap1 = 'snap1' ; snap3 = 'snap3' ; snap4 = 'snap4'	# snaparr formatted as [before,during,after]  for [short,medium,long] 
#									# accretion events. Uses dictionary to ensure that snaparr can be 										# manipulated into snaparr{i,j} array. DE05 incdices also listed in ###
d[snap0] = [131, 331]	### [200, 400] ###						
d[snap1] = [600, 1100, 2100] ### [669, 1169, 2169] ###
d[snap3] = [[100,150,200,250],[400,450,550,600],[1050,1150,1350,1450]]
	### [[00169,00219,00269,00319],[00469,00519,00619,00669],[01119,1219,01419,01519]] ###
d[snap4] = [[200,209,212,220],[860,875,880,900],[1525,1585,1595,1625]]
	### [[00269,00278,00281,00289],[00929,00944,00949,00969],[01594,01654,01664,01694]] ###
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - CONSTANTS AND CONVERSIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
pcAU = 206265.		# Conversion from parsec -> AU
AUm = 1.496e11		# Conversion from AU -> m
G = 6.67e-11		# Gravitational constant
Msol_kg = 1.998e30	# Conversion from Solar masses -> kg 
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
import os			# Standard Python modules 
import numpy as np
#
import r1_r			# Local modules
import sink_r
import pdisc_r
import rdisc_r
import pv_diag
import mass_comp
import exp_comp
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# Before beginning loops, remove old .dat files with radial exponents for SD and T
#
os.system('rm -r '+arch_dir+'SD_exps.dat '+arch_dir+'T_exps.dat')
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# For full analysis, loop over accretion tags as defined in array of L23 (ea_run)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# Loop over run number
#
for i in range(0, len(ea_run)):
#
	# Loop over ea inclination values
#
	for j in range(0, len(inclin)):
#
	# Loop over Keplerian restrictions
#
		for k in range(0, len(v_K)):
#
			print "\n Run "+str(ea_run[i])+ " from " \
			   "vK"+str(v_K[k])+"_"+str(inclin[j])+"i/ directory"+ \
			   " is now being analysed \n"
#
			plotdir = arch_dir+"PLOTS/"+str(ea_run[i])
			os.system('mkdir '+str(plotdir))
			plotdir = plotdir+"/"+"vK"+v_K[k]+"_"+inclin[j]+"i/"
			os.system('mkdir '+str(plotdir))
#
			arch_dir_tmp = arch_dir+"DATA/"+str(ea_run[i])+ \
			   "/"+"vK"+v_K[k]+"_"+inclin[j]+"i/"
#
			snaparr = d['snap'+str(ea_run[i])]
			snaparr = np.array(snaparr)
#
	# Read in rdisc.1 file
#
			hasharr_app, n_accr, r_d_kep, r_d_sigALMA, \
			   m_s, m_mri_d, m_d_kep, m_d_sigALMA, m_d_piv = \
			   r1_r.read(arch_dir_tmp, plotdir, ea_run[i], \
			   snaparr, v_K[k], inclin[j])
#
	# DE05.sink file(s) read (Planet sink data vs. radius vs. time)
#
			pmass, pradius, timearr = \
			   sink_r.read(arch_dir_tmp, plotdir, ea_run[i], snaparr, pcAU)
#
	# pdisc read (Disc parameters vs. radius for defined time)
#
			r, vkep = pdisc_r.read(arch_dir_tmp, plotdir, ea_run[i], \
			   snaparr, EA_timeref, EA_lenref, pmass, pradius, \
			   hasharr_app, n_accr, r_limit, r_d_kep, r_d_sigALMA,
			   timearr, v_K[k], inclin[j])
#
	# rdisc read (Disc parameters vs. time for defined radius)
#
			rdisc_r.read(arch_dir_tmp, plotdir, ea_run[i], hasharr_app, \
			   n_accr, r_inspec, v_K[k], inclin[j])
#
			if ((pv == "TRUE") and (inclin[j] != "0")):
#
	# Position-Velocity diagram plot
#
				pv_mass, kep_fit, raw_fit = pv_diag.pv(arch_dir_tmp, plotdir, ea_run[i], \
				   snaparr, v_K[k], inclin[j], r, vkep, EA_lenref, EA_timeref, \
				   pcAU, AUm, G, Msol_kg)	
#
	# Mass comparison of simulation to PV data
#
				if ( (raw_fit == "FALSE") and (kep_fit == "TRUE") and (mcomp_tmp == "TRUE")):
					mass_comp.comp(arch_dir_tmp, plotdir, ea_run[i], \
					   snaparr, timearr, v_K, inclin, m_s, m_mri_d, pmass, \
					   m_d_kep, m_d_piv, m_d_sigALMA, pv_mass, EA_lenref)
#
	# Convert all images to .eps format, suitable for publication
#
			os.system('for f in '+plotdir+'*.pdf; do pdftops -eps $f; done')
			os.system('rm -r '+plotdir+'*.pdf')
#
	# After looping over all selected ea runs, compare radial profile exponents for T and SD
#
if (exp == "TRUE"):
	plotdir = arch_dir + "/PLOTS/"
	exp_comp.comp(arch_dir, plotdir)
	os.system('for f in '+plotdir+'*.pdf; do pdftops -eps $f; done')
	os.system('rm -r '+plotdir+'*.pdf')
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Remove garbage from analysis directory and end program
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
os.system("rm -r *.pyc")
os.system("rm -r *~")
exit()
