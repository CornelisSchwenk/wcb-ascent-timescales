#---------------------------------------------------------
# Script to execute the WCB (Warm Conveyor Belt) selection.
# Created by Cornelis Schwenk on 12.12.2023
# Contact: c.schwenk@uni-mainz.de
#
# Note: This script should be executed on the MOGON-II supercomputer
#       after the merge_and_sum_NEST.slurm script has been run and
#       all trajectory files have been prepared.
#---------------------------------------------------------

# Include necessary scripts
include("WCB_selection.jl")

# Print start date and time for logging purposes
dt_start = Dates.now()
println("Start time: ", Dates.format(dt_start, "HH:MM:SS"))

# Write all meta files (WCB_tau_tstTIMESTEP_p00Y.nc)
# These files contain the index, ascent beginning time step, and t600 data.
write_all_WCB_meta_files()

# Re-include the selection script to load the files with glob("...")
include("WCB_selection.jl")

# Write WCB_tau_vars_realtime.nc
# This file contains all the data of the WCB trajectories along the entire simulation.
# No temporal re-ordering of the data is performed.
write_WCB_var_file_realtime()

#-------------------------------------------------
# This step must be done last!
# The function will include icon_WCB_functions.jl and
# requires WCB_tau_vars_realtime.nc to be available.
#-------------------------------------------------
# Write WCB_tau_vars_normed.nc
# This file contains the WCB ascent data stretched to a coherent ascent-time axis
# going from 0 (start) to 1 (end). This must be done after WCB_tau_vars_realtime.nc
# has been written.
write_tauWCB_var_file_normed()

# Print end date and time, and the duration of the script execution
dt_end = Dates.now()
println("End time: ", Dates.format(dt_end, "HH:MM:SS"))
println("Execution completed in $(round(dt_end - dt_start, Minute)) minutes.")
