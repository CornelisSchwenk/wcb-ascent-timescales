#!/bin/bash

#------------------------------------------------------------------------
#	This script is for merging and summating trajectory data output from ICON.
#	Created by Cornelis Schwenk 12.12.2023
#	c.schwenk@uni-mainz.de
#------------------------------------------------------------------------

#	!!NOTE!! You need CDO (climate data operator) to run this script

#------------------------------------------------------------------------
#	Explanation
#------------------------------------------------------------------------
#	Our Nested ICON trajectory data was so large, that we had to split it up
#	according to    1) variables split across 11 files (XX)
#			2) even then need multiple, up to four (Y)
#			3) domain (dom001, dom002 and dom003) (Z). 

#	The trajectory files thus look like this:
#		traj_XX_tstTIMESTEP_p00Y_dom00Z.nc

#	In the first step therefore have to MERGE all the traj_XX_* files.
#	For this we use CDO (Climate Data Operators from Max Planck Institute)
#	https://mpimet.mpg.de/cdo

#	In the second step we then SUMMATE over the domains. This is because
#	ICON outputs the data such that all trajectories are always written 
#	for each domain, but when a trajectory is for instance in domain 3, 
#	then the output for the domains 1 and 2 files will be zero at this
#	time step. So if a trajectory is in domain 1, then 2 and then 3 for
#	two time steps each, then the data for some variable var will look like

#		var_dom001 = [3.88,4.05,0.00,0.00,0.00,0.00,...]
#		var_dom002 = [0.00,0.00,3.98,4.01,0.00,0.00,...]
#		var_dom003 = [0.00,0.00,0.00,0.00,4.27,4.32,...]

#	and to obtain the final values we have to summate these arrays. This
#	is also done with CDO.

#	This script gets executed by a slurm caller in this folder: mogon_merge.slurm

#------------------------------------------------------------------------
#	Script
#------------------------------------------------------------------------
	#first enter the path to p_icon = path/to/icon/output/
	#and to p_traj = path/for/merged/trajectory/files/

p_icon="path/to/icon/output/"
p_traj="../nest_data/"

	#script goes to directory where icon data is stored
cd ${p_icon}

	#list all files that begin with traj_01. This will be one eleventh of all files.
	#We will loop over these and merge all 11 files.
	#The files will be called traj_total_tstTIMESTEP_p00Y_dom00Z.nc

	#--------------------------------
	#Merging
	#--------------------------------

l=$(ls traj_01*)

for i in $(echo $l)
do	
	name=$(echo $i | cut -c 9-31)	
	cdo merge ${p_icon}*$name* ${p_traj}traj_total_${name}.nc
done

	#--------------------------------
	#Summation
	#--------------------------------

	#Go to directory of merged trajectories

cd ${p_traj}

	#List all of the new files and then sum over all domains using CDO
	#The new files will be called traj_sum_tstTIMESTEP_p00Y.nc
l2=$(ls traj_total*dom001*)
for j in $(echo $l2)
do
	name=$(echo $j | cut -c 12-27)
	cdo enssum traj_total_${name}_dom001.nc traj_total_${name}_dom002.nc traj_total_${name}_dom003.nc traj_sum_${name}.nc
done

