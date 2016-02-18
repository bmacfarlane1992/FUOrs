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
# Date last modified: 17/02/2016 

arch_dir='/home/ben/Documents/WORK_PLANETS/PROJECTS/FUORS/'			# Achitecture directory

# Warning to ensure that correct file path entered into analyse_disc.f90

printf "\n Please ensure filepath correct for variables file in L153 of analyse_disc.f90 \n "
printf "\t Press [ENTER] to begin generation of data \n "
read ok

# Just in case, create DATA/ directory in project architecture

mkdir $arch_dir'DATA' 


# Loop over inclination values ($inclin$ within analyse_disc.f90 script) used in analyses

for i in 90 # 60 0
do

# Loop over Keplerian restrictions ($restrkep$ within analyse_disc.f90 script) used in analyses

	for j in 90 # 70
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

		for k in 3 # 4 
		do 

# Create directories for storage of output data files

			printf " \n EA run '$k' now being analysed \n \n "
			mkdir $arch_dir'DATA/'$k
			filedir=$arch_dir'DATA/'$k'/vK'$j'_'$i'i'
			mkdir $filedir
			mkdir $filedir'/pdisc/' $filedir'/cdisc/' $filedir'/rdisc/' $filedir'/pv_diag/'

# Now run executable over all DE05 files within relevant EA run directory.
 
			cd $arch_dir'ICs/'$k
			for x in DE05.du.0*
			do
				./../../GENERATE/seren/analysedisc "${x}"
			done

# Next, plot x-y distribution in SPLASH and send to appropriate directory

			cd $arch_dir'PLOTS/'
			dsplash ../DATA/$k'/vK'$j'_'$i'i/' -x 1 -y 2 -r 9 -dev SPLASH_xy.png
			mv SPLASH_xy.png $k'/vK'$j'_'$i'i/'

# Remove all garbage files and transfer useful data to pdisc/cdisc storage directory

			mv *.rdisc.1 $filedir'/rdisc/'
			rm -r *rdisc.1.*
			mv *.cdisc.1 $filedir'/cdisc/'
			mv *.pdisc.1 $filedir'/pdisc/'
			mv *hist_pv.1 *fit_Lpv.1 *fit_Rpv.1 $filedir'/pv_diag/'

		done
	done
done
