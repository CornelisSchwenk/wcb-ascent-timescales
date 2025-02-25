#!/bin/bash
#SBATCH -p smp
#SBATCH -A m2_jgu-tpchange
#SBATCH --time=20:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --constrain=skylake
#SBATCH --mem=100G
#SBATCH --job-name=calculate_wind_shear
#SBATCH --output=slurmfiles/%x_%j.out
#SBATCH --error=slurmfiles/%x_%j.err

module load lang/Julia
julia --project=. --threads 4 ../julia_code/calculate_mean_wind_shear.jl

#This will write onto the calculate_wind_shear.out file in the slurmfiles directory the arrays with the mean wind shear values in them. I copy-pasted those values into julia. I did not write them onto netCDF files because of laziness, values are correct nonetheless. 
