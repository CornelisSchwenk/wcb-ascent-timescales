#!/bin/bash
#SBATCH -p smp
#SBATCH -A m2_jgu-tpchange
#SBATCH --time=10:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH --job-name=Merge_and_sum_JOBNAME
#SBATCH --output=slurmfiles/%x_%j.out
#SBATCH --error=slurmfiles/%x_%j.err

#----------------------------------------------------------------
#	This script gets sent to mogon using the sbatch command.
#	In its current form, the merging (and summation) for the 
#	simulation as well as the entire selection and
#	WCB_tau_vars... file-writing will be done. This will take
#	take >20h.
#	Make sure to install all packages in the julia directory
#	first.
#----------------------------------------------------------------

#load CDO module
module load data/CDO/1.9.8-intel-2019b

#change mod to executable and execute for NEST
chmod +x merge_and_sum.sh
./merge_and_sum.sh 

#go to where the julia code is and execute the selection
#and WCB_tau_vars... file-writing process.
cd ../

#load julia
module load lang/Julia

cd ../julia_code/
julia --project=. --threads 4 do_WCB_selection.jl
julia --project=. write_PECR.jl 
