#
# pv_mass_dload.py
#
# Download batch of files for continuous plot comparing masses from PV diagram to simulation masses.
# Then run analyse_disc.f90 script and provide relevant pv fit files 
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
dload_dir = "/san/stellar3/dstamatellos/ea/"
file_res = 10
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
#
import os
import random
import math
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.axes as axes
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

def dload(arch_dir, plotdir, ea_run, snaparr, timearr, v_K, inclin, m_s, m_mri_d, pmass, m_d_kep, \
   m_d_piv, m_d_sigALMA, pv_mass, EA_lenref):
#
	# First, use timearr data to refine range of DE05 files to be downloaded, and generate string
	# to scp data
#
	# For EA [0, 1] runs
#
	if (snaparr.ndim == 1):
		file_n = len(snaparr)
		snaparr_tmp = [0]*file_n
		fcount = 0
		for i in range(0, len(snaparr)):
			snaparr_tmp[fcount] = snaparr[i] + 69
			fcount = fcount + 1
#
	# For EA [2, 3, 4, 5, 6] runs
#
	elif (snaparr.ndim == 2):
		file_n = len(snaparr)*len(snaparr[0])
		snaparr_tmp = [0]*file_n
		fcount = 0
		for i in range(0, len(snaparr)):
			for j in range(0, len(snaparr[0])):
				snaparr_tmp[fcount] = snaparr[i][j] + 69
				fcount = fcount + 1
#
	de05_low = snaparr_tmp[3] - 25
	de05_high = snaparr_tmp[8] + 25
	de05_n = (de05_high - de05_low) / file_res
#
	point_file = ""
	for i in range(0, de05_n):
		point_n = (i * file_res) + de05_low
		if (point_n < 1000):
			point_file = point_file+dload_dir+'DE05.00'+str(point_n)+' '
		if (point_n > 1000):
			point_file = point_file+dload_dir+'DE05.0'+str(point_n)+' '
#
	# Now download
#
	os.system("scp bmacfarlane@stargate.uclan.ac.uk:'"+point_file+"' "+arch_dir+"ICs/"+ea_run"/pv_cont/")
