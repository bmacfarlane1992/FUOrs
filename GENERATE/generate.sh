#!/bin/sh
#
# Script to run analysis of episodic accretion runs based on parameters of disc definition vs.
# Keplerian profile and inclination. Parameters defined within ../seren/src/analyse/analyse_disc.f90
#
# Script outputs relevant pdisc and cdisc files into appropriate directories
# for further analysis in FUOrs project
#
# Author: Ben MacFarlane
# Contact: bmacfarlane@uclan.ac.uk
# Date last modified: 29/01/2016 

arch_dir='/home/ben/Documents/WORK_PLANETS/PROJECTS/FUORS/'			# Achitecture directory

# Warning to ensure that correct file path entered into analyse_disc.f90

printf "\n Please ensure filepath correct for variables file in L153 of analyse_disc.f90 \n "
printf "\t Press [ENTER] to begin generation of data \n "
read ok

# Loop over inclination values ($inclin$ within analyse_disc.f90 script) used in analyses

for i in 90 75 60 30 0
do

# Loop over Keplerian restrictions ($restrkep$ within analyse_disc.f90 script) used in analyses

	for j in 90 70 50
	do

# Write variables to file for reading within analyse_disc.f90

		rm -r $arch_dir'GENERATE/varfile.dat'
		varfile=$arch_dir'GENERATE/varfile.dat'
		echo "$i" >> "$varfile"
		echo "$j" >> "$varfile"

# Remake seren analysedisc script

		cd $arch_dir'GENERATE/seren/'
		make clean
		make analysedisc

# Loop over EA runs to be analysed

		for k in 3 4 5 6
		do 

# Create directories for storage of output data files

			printf " \n EA run '$k' now being analysed \n \n "
			filedir=$arch_dir'DATA/'$k'/vK'$j'_'$i'i'
			rm -r $filedir'/pv_diag/'*raw*
			mkdir $filedir
			mkdir $filedir'/pdisc/' $filedir'/cdisc/' $filedir'/rdisc/' $filedir'/pv_diag/'

# Now run executable over all DE05 files within relevant EA run directory.
 
			cd $arch_dir'ICs/'$k
			for x in DE05.du.0*
			do
				./../../GENERATE/seren/analysedisc "${x}"
			done

# Remove all garbage files and transfer useful data to pdisc/cdisc storage directory

			mv *.rdisc.1 $filedir'/rdisc/'
			rm -r *rdisc.1.*
			mv *.cdisc.1 $filedir'/cdisc/'
			mv *.pdisc.1 $filedir'/pdisc/'
			mv *.raw_pv.1 *hist_pv.1 *fit_Lpv.1 *fit_Rpv.1 $filedir'/pv_diag/'

		done
	done
done
