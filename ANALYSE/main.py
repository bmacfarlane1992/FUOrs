#
# main.py
#
# Python program to execute read and analysis of simulation files.
#
# Author: Benjamin MacFarlane
# Date: 30/01/2016
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
v_K = ["90"] # ["90","70","50"]						# Keplerian velocity restriction on disc mass/radius
inclin = ["90"] # ["0","30","60","75","90"]				# Inclination of disc being analysed
#
batch = "TRUE"							# Choose whether ("TRUE") or not ("FALSE") to batch
batch_accrun = [3] # [3,4,5,6]					# process runs
#
acc_run = 3						# If batch = "FALSE", carry out analysis of single run (debugging)			
#
r_limit = 150							# Limit of radial plots in pdisc analyses
spline = "TRUE" 						# Choose whether or not to smooth surface density distributions
#
r_inspec = 150 							# Radius at which disc parameters inspected in rdisc
#
pv_plot = "TRUE"						# Choose whether ("TRUE") or not ("FALSE") to generate PV diagram
#
EA_timeref = ["BEFORE","DURING","AFTER"]			# Define names of EA snapshots for EA length reference	
EA_lenref = ["SHORT","MEDIUM","LONG"]				# and time reference to EA outburst event
#
d = {}
snap0 = 'snap0' ; snap1 = 'snap1' ; snap2 = 'snap2'		# snaparr formatted as [before,during,after]  for [short,medium,long] 
snap3 = 'snap3' ; snap4 = 'snap4' ; snap5 = 'snap5' ; snap6 = 'snap6'	# accretion events. Must use dictionary to ensure that snaparr
#									#  can be manipulated into snaparr{i,j} array
d[snap0] = [100, 500]						
#		DE05 ext for ea 0: [169, 569] 			# Respective DE05 indices also listed
d[snap1] = [100, 1100, 2100]
#		DE05 ext for ea 1: [169, 1169, 2169]
d[snap2] = [[150,190,230],[800,840,880],[1500,1545,1590]]
#		DE05 ext for ea 2: [[00219,00259,00299],[00869,00909,00949],[01569,01614,01659]]
d[snap3] = [[100,175,250],[400,500,600],[1050,1250,1450]]
#		DE05 ext for ea 3: [[00169,00244,00319],[00469,00569,00669],[01119,01319,01519]] 
d[snap4] = [[200,210,220],[860,880,900],[1525,1575,1625]]
#		DE05 ext for ea 4: [[00269,00279,00289],[00929,00949,00969],[01594,01644,01694]]
d[snap5] = [[170,190,210],[740,760,780],[1360,1390,1420]]
#		DE05 ext for ea 5: [[00239,00259,00279],[00809,00829,00849],[01429,01459,01489]] 
d[snap6] = [[185,200,215],[840,870,900],[1560,1600,1640]]
#		DE05 ext for ea 6: [[00254,00269,00284],[00909,00939,00969],[01629,01669,01709]]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
	# Import python library modules
#
import os
import numpy as np
#
	# Import local modules
#
import r1_r
import sink_r
import pdisc_r
import rdisc_r
import pv_diag
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
# Loop over inclination values
#
for i in range(0, len(inclin)):
#
# Loop over Keplerian restrictions
#
	for j in range(0, len(v_K)):
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# For full analysis, loop over accretion tags as defined in array of L23 (batch_accrun) #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
		if (batch == "TRUE"):
#
			for k in range(0, len(batch_accrun)):
				print "\n Run "+str(batch_accrun[k])+ " from " \
				   "vK"+str(v_K[j])+"_"+str(inclin[i])+"i/ directory"+ \
				   " is now being analysed \n"
#
				plotdir = arch_dir+"PLOTS/"+str(batch_accrun[k])
				os.system('mkdir '+str(plotdir))
				plotdir = plotdir+"/"+"vK"+v_K[j]+"_"+inclin[i]+"i/"
				os.system('mkdir '+str(plotdir))
#
				batch_archdir = arch_dir+"DATA/"+str(batch_accrun[k])+ \
				   "/"+"vK"+v_K[j]+"_"+inclin[i]+"i/"
#
				snaparr = d['snap'+str(batch_accrun[k])]
				snaparr = np.array(snaparr)
#
	# Read in rdisc.1 file
#
				acc_count, acctag_num, hasharr_app, n_accr, r_d_kep, \
				   m_s, m_d_kep, m_d_sigALMA, m_d_piv = \
				   r1_r.read(batch_archdir, plotdir, batch_accrun[k], \
				   snaparr, v_K[j], inclin[i])
#
	# DE05.sink file(s) read (Planet sink data vs. radius vs. time)
#
				pmass, pradius, timearr = \
				   sink_r.read(batch_archdir, plotdir, batch_accrun[k], snaparr)
#
	# pdisc read (Disc parameters vs. radius for defined time)
#
				pdisc_r.read(batch_archdir, plotdir, batch_accrun[k], \
				   snaparr, EA_timeref, EA_lenref, pmass, pradius, \
				   hasharr_app, n_accr, r_limit, spline, r_d_kep, timearr, v_K[j])
#
	# rdisc read (Disc parameters vs. time for defined radius)
#
				rdisc_r.read(batch_archdir, plotdir, batch_accrun[k], hasharr_app, \
				   n_accr, r_inspec, v_K[j], inclin[i])
#
	# Position-Velocity diagram plot
#
				if (pv_plot == "TRUE"):
					pv_diag.pv(batch_archdir, plotdir, batch_accrun[k], \
					   snaparr, v_K[j], inclin[i], EA_lenref, EA_timeref, \
					   m_s, m_d_kep, m_d_piv, m_d_sigALMA)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# De-bugging analysis - requires entry of both acc_num (L19) and run_tag (L20) #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
		else:
#
			plotdir = arch_dir+"PLOTS/"+str(acc_run)
			os.system('mkdir '+str(plotdir))
			plotdir = plotdir+"/"+"vK"+v_K[j]+"_"+inclin[i]+"i/"
			os.system('mkdir '+str(plotdir))
#
			arch_dir = arch_dir+"DATA/"+str(acc_run)+"/"+"vK"+v_K[j]+"_"+inclin[i]+"i/"
#
			snaparr = d['snap'+str(acc_run)] ; snaparr = np.array(snaparr)
			snaparr = np.array(snaparr)
			print "\n Run "+str(acc_run)+"from vK"+str(v_K[j])+"_"+str(inclin[i])+ \
			   "i/ directory is now being analysed \n"
#
			acc_count, acctag_num, hasharr_app, n_accr, r_d_kep, \
			   m_s, m_d_kep, m_d_sigALMA, m_d_piv = \
			   r1_r.read(arch_dir, plotdir, acc_run, snaparr, v_K[j], inclin[i])
#
			pmass, pradius, timearr = \
			   sink_r.read(arch_dir, plotdir, acc_run, snaparr)
#
			pdisc_r.read(arch_dir, plotdir, acc_run, snaparr, EA_timeref, \
			   EA_lenref, pmass, pradius, hasharr_app, n_accr, \
			   r_limit, spline, r_d_kep, timearr, v_K[j])
#
			rdisc_r.read(arch_dir, plotdir, acc_run, hasharr_app, n_accr, r_inspec, \
			   v_K[j], inclin[i])
#
			if (pv_plot == "TRUE"):
				pv_diag.pv(arch_dir, plotdir, acc_run, snaparr, v_K[j], inclin[i], \
				   EA_lenref, EA_timeref, m_s, m_d_kep, m_d_piv, m_d_sigALMA)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# Remove garbage from analysis directory and end program #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
os.system("rm -r *.pyc")
os.system("rm -r *~")
exit()
#
#
