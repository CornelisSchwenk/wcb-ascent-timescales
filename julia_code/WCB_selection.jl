#---------------------------------------------------------------------------
#       "Skript for the selection of WCB trajectories from ICON trajectory files"
#       "Created by Cornelis Schwenk 22.11.2023"
#       "c.schwenk@uni-mainz.de"
#---------------------------------------------------------------------------

#       ----------   INSTRUCTIONS   ----------

#       "to write the metdata files just include this script in julia and      " #
#       "then exectute the function: write_all_WCB_meta_files()              " #

#       "to write the large WCB files with variables just include this script  " #
#       "in julia and then exectute the function: write_WCB_var_file()       " #            

#       -------------   REMARKS --------------

#	"The Eulerian data from ICON (the icon output) is not included in the"
#	"Zenodo dataset. To get this contact me at the email above"

#--------------------------------------------------------

include("use_config.jl")
using Glob, NCDatasets, Plots
using Plots.PlotMeasures
using Statistics
using BenchmarkTools, Dates, DataStructures
using DelimitedFiles
using NaNMath; nm=NaNMath
using Interpolations
using Base.Threads

#--------------------------------------------------------
#       "define paths to ICON data			"
#--------------------------------------------------------

# "path to trajectory from nested run"
global nest_path = "../nest_data/"			

# "path of ICON Eulerian data. not included in Zenodo dataset"
global nest_NWP_path = "path/to/icon/output/"	

# "list of all trajectory variables"
varlist = ["idx", "CellId", "lon", "lat", "alt", "rtime", "w_v", "rho", "t", "p", "pv", "tmin", "tmax","qv", "qc", "qr", "qi", "qs", "qg", "qh", "qnc", "qnr", "qni", "qns", "qng", "qnh", "qi_in", "qs_in", "qg_in", "qh_in", "qr_in", "qi_out", "qs_out", "qg_out", "qh_out", "qr_out","r_evap", "f_evap", "evap", "subl", "r_melt", "c_melt","cond", "depo", "qc_nuc", "qi_hom", "qihh", "wbf", "qx_dep", "satad2", "r_frez", "qx_rim","ni_hom", "ni_hh", "nc_nuc", "n_sim","qvturb", "qcturb", "qiturb", "qvconv", "qcconv", "qrconv", "qiconv", "qsconv","qvturc","qcturc"]

# "list of microphysical traj. variables"
# "separate lists in case different output"
microlist = ["qi_in", "qs_in", "qg_in", "qh_in", "qr_in", "qi_out", "qs_out", "qg_out", "qh_out", "qr_out", "r_evap", "f_evap", "evap", "subl", "r_melt", "c_melt", "cond", "depo", "qc_nuc", "qi_hom", "qihh", "wbf", "qx_dep", "satad2", "r_frez", "qx_rim", "ni_hom", "ni_hh", "nc_nuc", "n_sim", "qvturb", "qcturb", "qiturb","qvconv", "qcconv", "qrconv", "qiconv", "qsconv","qvturc","qcturc"]

#-----------------------------------
#       "ICON trajecotry files."
#	"This list is empty if ICON data is not downloaded and the script merge_and_sum.sh from"
#	"EXECUTE_WCB_SELECTION is not executed"
#----------------------------------

nest_fls = glob(nest_path*"traj_sum*")[1:end-12] #"the last 12 don't have enough hours left in sim"
nest_WCB_fls = glob(nest_path*"WCB_tau_tst*")

#-----------------------------------------------------------------------
#       "define important indices etc."
#-----------------------------------------------------------------------

	#"!! make nest_skp and glob_skp such that time resolution is 30 minutes !!"
global nest_start = 1			#"where to start in traj files, nest"
global nest_skp = 4			#"how many data points to skip"

#"nested trajectory files must be summated (done in different skript), because"
#"values in _domX.nc are zero if traj not in domain X. So one must summate"
#"all trajectory files from different domains. The time however must then be"
#"divided by the amount of domains, because it is always written, hence nest_factor = 3"
global nest_offset_nan = 0		#"if NaNs are -999, then offset must be 999"
global nest_factor     = 3		


#"for writing the names in the netcdf file"
#"This is purely for name writing and if dependent on how trajectory files are named!"
global nest_name_ind = 18
global nest_title_ind = 27

#----------------------------------------------------------------------
#       "THE THREE FOLLOWING FUNCTIONS ARE THE IMPORTANT ONES           " #
#       "THEY USE ALL OF THE OTHER FUNCTIONS DEFINED BELOW               " #
#----------------------------------------------------------------------

function write_all_WCB_meta_files()
        #"THIS IS A BIG FUNCTION. IT WRITES ALL OF THE METADATA WCB FILES."
	#"Meta WCB files give the index and positions of wcbs within the trajectory dataset."
        
	#"This function goes through each trajectory file and writes a separate ncfile that"
	#"contains the information on:" 
	#"	- index of trajs that are WCB (w)"
	#"	- where they begin their ascent (a)" 
	#"	- how long the ascent takes in time steps (t6)."
	#"This is done by the "write_WCB_ncfile" function, see for more detail on selection process"

        fl_list = nest_fls
	println("length fl_list = $(length(fl_list))")
        for i in 1:length(fl_list) 	
	    f = fl_list[i]
	    print("$(i), ")
            write_WCB_ncfile(f)
        end
	println("! done !")
end

function write_WCB_var_file_realtime()
	#"This function does almost the same as write_WCB_var_file, except rtime is not rearanged"
	#"and the data is written all the way to the end of the simulation."
	#"All trajs begin at first runtime output. If trajs begin later (trajs might be started"
	#"at diffferent rtimes), then they are NaN until data is written. Example:"
	#"	traj 901 starts at rtime=0, then some_var = [1,2,3,4,...,192]"
	#" 	if traj 4000 starts at rtime=4, then some_var = [NaN,NaN,NaN,NaN,5,6,...,192]"

	#"The file WCB_tau_vars_realtime.nc is written."
	#"This uses primarily the function get_var_realtime"

        pth = nest_path

	ntrajs = length(get_w_a_all()[1])
        ds_out_name = "WCB_tau_vars_realtime.nc"
        ds_out = Dataset(pth*ds_out_name,"c")
        defDim(ds_out,"time",192)
        defDim(ds_out,"ntraj",ntrajs)
        ds_out.attrib["title"] = "WCB variables in real time"

        vardict = Dict()
        ncdict = Dict()
	println("------------STARTING------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
        for i in varlist[varlist .!= "rtime"]
                println("writing variable $(i)")
                ncdict[i] = defVar(ds_out,i,Float64,("ntraj","time"))
                ncdict[i][:,:] = get_var_realtime(i)
		GC.gc()
        end
	
	w,a,t6 = get_w_a_all()
	t6_ = defVar(ds_out,"tau_600",Float64,("ntraj",))	
	t6_.attrib["units"] = "hours"
	t6_[:] = t6 ./ 2
	
	r_t = defVar(ds_out,"rtime",Float64,("time",))
	r_t[:] = collect(0.5:0.5:96)

	dp2 = find_max_dp2()
	dp2_ = defVar(ds_out,"dp2",Float64,("ntraj",))
	dp2_.attrib["units"] = "hPa"
	dp2_[:] = dp2

        close(ds_out)
	println("------------FINISHED------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
end

function write_tauWCB_var_file_normed()
        #" !! This function must be executed AFTER write_WCB_var_file_realtime !!"
	#"This function won't run is WCB_tau_vars_realtime.nc doesn't exist already"

	#"This function creates the WCB_tau_vars_normed.nc file." 
	#"It writes variables in a normalized ascent time scale from 0 to 1 in 101 steps."
	#"0 is beginning of ascent. 1 is end of ascent. It takes the data from WCB_vars_realtime but"
	#"stretches/contract the ascent part we want such that all arrays have the same length."
	#"The times are still such that the first time entry is at the beginning of the ascent."
	#"rtime is still written."
	#"There are no more NaNs. If we take the example from the function above:"
	#" 	traj 901 has t6 = 2h. Then some_var = [1,2,3,4,NaN,NaN,...,NaN]"
	#" 			       and some_var_normed = [1,1.2,1.4,...,3.8,4]"

	#"the function "get_var_tWCB_normed" is the one that is primarily used here"
	
	pth = nest_path
	
	if length(glob(pth*"WCB_tau_vars_realtime.nc")) == 0
		error("You must create the WCB_tau_vars_realtime.nc file first, or
		      download it from zenodo and place it in ../nest_data/")
	end
	
	include("icon_WCB_functions.jl")

	global t6n = dsnr["tau_600"][:,:]
	t6 = t6n
	global pnr = dsnr["p"][:,:] ./ 100 
	starts = [find_start(k) for k in 1:length(t6)]
        stops = [find_stop(k) for k in 1:length(t6)]

	ntrajs = length(t6)
        ds_out_name = "WCB_tau_vars_normed.nc"
        ds_out = Dataset(pth*ds_out_name,"c")
        defDim(ds_out,"time",101)
        defDim(ds_out,"ntraj",ntrajs)
        ds_out.attrib["title"] = "variables during normed WCB ascent time"

        vardict = Dict()
        ncdict = Dict()
	println("------------STARTING------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
	for i in varlist
		print("$(i), ")
		println("writing variable $(i)")
		ncdict[i] = defVar(ds_out,i,Float64,("ntraj","time"))
		ncdict[i][:,:] = get_var_tWCB_normed(i,starts,stops)
	end
	t6_ = defVar(ds_out,"tau_600",Float64,("ntraj",))	
	t6_.attrib["units"] = "hours"
	t6_[:] = t6

	tWCB = stops .- starts
	tWCB_ = defVar(ds_out,"tau_WCB",Float64,("ntraj",))	
	tWCB_.attrib["units"] = "hours"
	tWCB_[:] = tWCB ./ 2
        
	dp2 = dsnr["dp2"][:]
	dp2_ = defVar(ds_out,"dp2",Float64,("ntraj",))
	dp2_.attrib["units"] = "hPa"
	dp2_[:] = dp2
	
	close(ds_out)
	println("------------FINISHED------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
end
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#"       THE FOLLOWING FUNCTIONS ARE THE ONES USED BY THE FOUR BIG FUNCTIONS ABOVE    "#
#----------------------------------------------------------------------

function write_WCB_ncfile(fl::String)
        #"iterate this function over the file list of the trajectories"
        #"for f in nfl write_WCB_ncfile(f) end"

        #"this write only the small files with the metadata of where to find the WCB trajs"
        #"within the traj files"

	#"The selection is done by the funciton "WCB_sel" that returns index (w), ascent start (a)"
	#"and ascent time (t600). See it for more details on selection."

	start = nest_start
	if fl[end-18:end-8] == "tst00000003"
		start = 4
	elseif fl[end-18:end-8] == "tst00001439"
		start = 5
	else
		start = 1
	end

	nan_offset = nest_offset_nan
	files = nest_fls[1:length(nest_WCB_fls)]
	name_ind = nest_name_ind
	title_ind = nest_title_ind
	pth = nest_path
	factor = nest_factor
	skp = nest_skp
	
        println("going through file $(fl)")
        ds = Dataset(fl)
        wcb,asc,tau = WCB_sel(ds,start,nan_offset,factor,skp)
	
	#"If file contains any WCB trajectories, write netCDF file."
        if length(wcb) > 0
                ds_out_name = fl[end-name_ind:end]
                ds_out = Dataset(pth*"WCB_tau_$(ds_out_name)","c")
                defDim(ds_out,"pos",size(wcb)[1])
                ds_out.attrib["title"] = "WCB positions for $(fl[end-title_ind:end])"

                pos = defVar(ds_out,"position",Int64,("pos",), attrib = OrderedDict("units" => "Row number in file"))
                ascent = defVar(ds_out,"ascent",Int64,("pos",), attrib = OrderedDict("units" => "Array pos. at which WCB ascent begins"))
                tau600 = defVar(ds_out,"tau600",Int64,("pos",), attrib = OrderedDict("units" => "Time steps until completion of ascent"))
				pos[:] = wcb
        		ascent[:] = asc
				tau600[:] = tau
                println("found $(size(wcb)) WCBs")
                close(ds_out)
        elseif length(wcb) == 0
                println("found no WCB trajectories, not writing file")
        end
	close(ds)
end

function find_tau_and_asc(p)
	#"This function takes a 1 dim pressure array in hPa and returns ascent time"
	#"and where ascent begins (so asc and tau600)."
        #"p must be in hPa and a 1 dim array"
        DT = Inf
        t = []
        for i in 1:(length(p)-1)
		for j in i:length(p)
		       dp = p[i] - p[j]
		       dt = j-i
		       if dp > 600
			       if dt < DT
				       push!(t,(i,dt))
				       DT=dt
			       end
			       break
		       end
	       end
        end
        return t[end]
end


function find_tau_x00(p,x)
	#"This function takes a 1 dim pressure array in hPa and returns ascent time"
        #"p must be in hPa and a 1 dim array"
        DT = Inf
        t = []
        for i in 1:(length(p)-1)
		for j in i:length(p)
		       dp = p[i] - p[j]
		       dt = j-i
		       if dp > x
			       if dt < DT
				       push!(t,dt)
				       DT=dt
			       end
			       break
		       end
	       end
        end
        return t[end]
end

function find_WCB(ds,start,nan_offset,factor,skp)	
	#"This function finds WCBs in a netCDF dataset of trajectories"

	#"convert rtime into hours"
	rt = ds["rtime"][start:skp:end] ./ (3600 * factor)
	#"convert pres into hPa"
        pres = (ds["p"][1:1:end,start:skp:end] .+ nan_offset)./ 100

        wcb_p = Int[]
	asc = Int[]
	tau = Int[]
        
	#"find all trajectories that go up 600hPa in at least 48h"
        for i in 1:(size(rt)[1]-Int(48*2))
                p0 = pres[:,i]
                p48 = pres[:,i+Int(48*2)]
		
                dp = p0 .- p48
                wcbs = findall(dp .> 600)
                for x in wcbs
                        push!(wcb_p, x)
                end
        end

	#"sort arrays and make sure no double couting"
        uniq = findfirst.(isequal.(unique(wcb_p)),[wcb_p])
        wcb_sort = sortperm(wcb_p[uniq])

	wcb_out = wcb_p[uniq][wcb_sort]

	for i in 1:length(wcb_out)
		a,t600 = find_tau_and_asc(pres[wcb_out[i],:])
		push!(asc,a)
		push!(tau,t600)
	end

	return (wcb_out,asc,tau)
end

function WCB_sel(ds,start,nan_offset,factor,skp)
	#"This function takes the WCB ascents from "find_WCB" function and then"
	#"does a second selection process to make sure they are part of the cyclone"

        println("---------------------------------------------------")
        println("""going through dataset for $(ds.attrib["time"])""")

        rt = ds["rtime"][start:skp:end] ./ (3600 * factor)

	#"This first rtime value behaves weirdly, so use second one minus half hour"
        rt0 = rt[2] .- 0.5
        println("there will be $(size(rt)[1]) loops - $(rt[2] .- 0.5) to $(rt[end])")

        wcbs,asc,tau = find_WCB(ds,start,nan_offset,factor,skp)

        lat = (Array(ds["lat"])[wcbs,start:skp:end] .+ nan_offset) .* 180/pi
        lon = (Array(ds["lon"])[wcbs,start:skp:end] .+ nan_offset) .* 180/pi

        cyc_mask = Int.(ones(size(lat)[1]))

	#----------------------------------
	#"The second selection is a visual confirmation for specific time." 
	#"In our case we use 20170922_08 and define a boundary."
	#"If trajs outside of boundary at this time, mask is set to zero (or false)."
	#----------------------------------
	t_check = 56 #"rtime of 20170922_08 is 56"
	j_check = Int((t_check - (rt[2] - 0.5))*2) + 1
	t_check2 = 80 #"rtime of 20170923_08 is 80"
	j_check2 = Int((t_check2 - (rt[2] - 0.5))*2) + 1
	temp_cyc_mask2 = Int.(ones(size(lat)[1]))
	for i in 1:size(lat)[1]
		if !(-52 < lon[i,j_check] < -5) .| !(33 < lat[i,j_check] < 60)
			temp_cyc_mask2[i] = 0
		elseif (lon[i,j_check2] < -25) .& (lat[i,j_check2] .< 40)
			temp_cyc_mask2[i] = 0
		end
	end
	cyc_mask = cyc_mask .& temp_cyc_mask2

        println("")
        println("---------------------------------------------------")
        return (wcbs[Bool.(cyc_mask)],asc[Bool.(cyc_mask)],tau[Bool.(cyc_mask)])
end

function WCB_info(fl::String)
        #"This function gives you the WCB metadata for the trajectory file"
    	name_ind = nest_name_ind
    	pth = nest_path

        nm = fl[end-name_ind:end]
        wcbfl = pth*"WCB_tau_$(nm)"
		ds = Dataset(wcbfl)
        return (Array(ds["position"]),Array(ds["ascent"]),Array(ds["tau600"]))
	close(ds)
end

function get_w_a_all()
	#"This function returns the WCB metadata for ALL trajectories."

	files = [glob("traj_sum*$(nf[end-15:end])",nest_path)[1] for nf in nest_WCB_fls]

	w_tot = Vector{Float64}()
	a_tot = Vector{Float64}()
	t6_tot = Vector{Float64}()
	for f in files
		w,a,t6 = WCB_info(f)
		w_tot = vcat(w_tot,w)
		a_tot = vcat(a_tot,a)
		t6_tot = vcat(t6_tot,t6)
	end
	return Int.(w_tot),Int.(a_tot),Int.(t6_tot)
end


function find_start(k)
	#"This function finds where the WCB ascent starts before t600 begins."
	#"It looks for where the pressure velocity is on average smaller than"
	#"8hPa per hour in the hour before t600"
	if (@isdefined pnr)
		pk = pnr[k,:]
		a_s = an
		t_s = tn
		n_s = ns
	elseif (!@isdefined pnr)
		error("pnr must be defined globally first")
	end

	start_max = find_start_max(k,pk,a_s,n_s,t_s)
	start_min = find_start_min(k,pk,a_s,n_s,t_s)
	
	if start_max > start_min
		startlist = start_min:start_max
	elseif start_max < start_min
		startlist = start_max:start_min
	elseif start_max == start_min
		startlist = start_max
	end
	
	pstartlist = [pk[i] for i in startlist]
	startarg = argmax(pstartlist)
	return startlist[startarg]
	
end

function find_start_min(k,pk,a_s,n_s,t_s)
	#"This function finds the minimum time after which WCB ascent no"
	#"longer is larger than 8hPa/h after t600 is completed."
	#"It goes backwards from the end and looks if the mean ascent is"
	#"faster than 8hPa per hour. This function is used by find_stop"
	a_start = a_s[k] + n_s[k]
	if a_start == 1
		return 1
	end
	for i in 1:a_start
		p_diff_min = -4 .* (a_start - i)
		if i == a_start
			return a_start
		elseif pk[a_start] .- pk[i] .< p_diff_min
			return i
		end	
	end
end

function find_start_max(k,pk,a_s,n_s,t_s)
	#"This function finds where the WCB ascent starts before t600 begins."
	#"It looks for where the pressure velocity is on average smaller than"
	#"8hPa per hour in the hour before t600"
	if a_s[k] == 1
		return a_s[k] + n_s[k]
	elseif a_s[k] == 2
		if (pk[n_s[k] + a_s[k]] - pk[n_s[k] + a_s[k] - 1]) <= -4
			return n_s[k] + a_s[k] - 1
		else
			return n_s[k] + a_s[k]
		end
	else
		counter=0
		for i in reverse(n_s[k]+3:n_s[k]+a_s[k])
			counter = i-2
			if (pk[i] - pk[i-2]) .> -8
				if (pk[i] - pk[i-2]) .<= 0
					return i-2
				elseif (pk[i] - pk[i-1]) .<= 0
					return i-1
				else
					return i
				end
			end
		end
		if counter == n_s[k] + 1
			return n_s[k] + 1
		end
	end
end

function find_stop(k)
	#"This function finds where the WCB ascent stops after t600 is completed."
        #"It uses find_stop_max and find_stop_min and looks for where in between"
	#"those two the pressure is minimal"
	if (@isdefined pnr)
		pk = pnr[k,:]
		a_s = an
		t_s = tn
		n_s = ns
	elseif (!@isdefined pnr)
		error("pnr must be defined globally first")
	end
	
	stop_max_1 = find_stop_max_1(k,pk,a_s,n_s,t_s)
	stop_max_2 = find_stop_max_2(k,pk,a_s,n_s,t_s)
	stop_max = minimum([stop_max_1,stop_max_2])	
	stop_min = find_stop_min(k,pk,a_s,n_s,t_s)

	if stop_max > stop_min
		stoplist = stop_min:stop_max
	elseif stop_max < stop_min
		stoplist = stop_max:stop_min
	elseif stop_max == stop_min
		stoplist = [stop_max]
	end
	pstoplist = [pk[i] for i in stoplist]
	stoparg = argmin(pstoplist)
	return stoplist[stoparg]

end

function find_stop_max_1(k,pk,a_s,n_s,t_s)
	#"This function finds the minimum time after which WCB ascent no"
	#"longer is larger than 8hPa/h after t600 is completed."
	#"It goes backwards from the end and looks if the mean ascent is"
	#"faster than 8hPa per hour. This function is used by find_stop"
	t6stop = a_s[k] + t_s[k] + n_s[k]
	k_to_stop = 192 - t6stop
	if k_to_stop == 0
		return 192
	end
	for i in k_to_stop:-1:0
		p_diff_min = -4 .* i
		if i == 0
			return t6stop
		elseif pk[t6stop + i] .- pk[t6stop] .< p_diff_min
			return t6stop + i
		end	
	end
end

function find_stop_max_2(k,pk,a_s,n_s,t_s)
	#"This function finds the minimum time after which WCB ascent no"
	#"longer is larger than 8hPa/h after tWCB is completed."
	#"It goes backwards from the end and looks if the mean ascent is"
	#"faster than 8hPa per hour. This function is used by find_stop"
	tWCBstop = find_stop_min(k,pk,a_s,n_s,t_s)
	k_to_stop = 192 - tWCBstop
	if k_to_stop == 0
		return 192
	end
	for i in k_to_stop:-1:0
		p_diff_min = -4.5 .* i
		if i == 0
			return tWCBstop
		elseif pk[tWCBstop + i] .- pk[tWCBstop] .< p_diff_min
			return tWCBstop + i
		end	
	end
end

function find_stop_min(k,pk,a_s,t_s,n_s)
	#"This function finds the minimum time after which WCB ascent no"
	#"longer is larger than 8hPa/h after t600 is completed."
	#"It looks for where the pressure velocity is on average smaller than"
	#"8hPa per hour in the hour after t600. This function is used by find_stop."

	if a_s[k] + t_s[k] + n_s[k] == 190
		if (pk[192] - pk[190]) .<= -8
			return 192
		elseif (pk[191] - pk[190]) .<= -4
			return 191
		else
			return 190
		end
	elseif a_s[k] + t_s[k] + n_s[k] == 191
		if (pk[192] - pk[191]) .<= -4
			return 192
		else
			return 191
		end
	elseif a_s[k] + t_s[k] + n_s[k] == 192
		return 192
	else
		counter = 0
		for i in (a_s[k]+t_s[k]+n_s[k]):190
			counter=i+2
			if (pk[i+2] - pk[i]) .> -8
				if (pk[i+2] - pk[i]) .<= 0
					return i+2
				elseif (pk[i+1] - pk[i]) .<= 0
					return i+1
				else
					return i
				end
			end
		end
		if counter == 192
			return counter
		end
	end
end

function get_var_tWCB_normed(vr::String,starts,stops)

	#" !! This function can only be executed once WCB_tau_vars_realtime.nc has been written !!" 

	fl = nest_path*"WCB_tau_vars_realtime.nc"
	if (!@isdefined t6n)
		error("t6n must be defined globally first!") 
	elseif (@isdefined t6n)
		t6 = t6n
	end
	
	#"Initialize array"
        v_tot = Vector{Matrix{Float64}}()

	#"Load data from WCB_tau_vars.nc"
	ds = Dataset(fl)
	vr_ar = Array(ds[vr])

	v_f = fill(NaN,(length(t6),101))
	
	#"Loop over all trajectories. LinearInterpolation is used to remap all variable array"
	#"from WCB_tau_vars_realtime.nc onto the size 101."
	if vr != "rtime"
		for i in 1:length(t6)
			v = vr_ar[i,starts[i]:stops[i]]
			f_int = LinearInterpolation(1:length(v),v,extrapolation_bc=Line())
			stps = (length(v)-1)/100
			upto = length(v)
			v_f[i,1:end] = f_int(1:stps:upto)
		end
	elseif vr == "rtime"
		for i in 1:length(t6)
			v = vr_ar[starts[i]:stops[i]]
			f_int = LinearInterpolation(1:length(v),v,extrapolation_bc=Line())
			stps = (length(v)-1)/100
			upto = length(v)
			v_f[i,1:end] = f_int(1:stps:upto)
		end
	end
	v_tot = vcat(v_tot,[v_f[i,:] for i in 1:size(v_f)[1]])	
	close(ds)
        v_out = reduce(vcat,transpose.(v_tot))
	if vr in microlist
		v_out = v_out .- v_out[:,1]
	end
	return v_out
end

function get_var_realtime(vr::String)
	#"This function retrieves a variable for all WCB trajectories by going through all"
	#"files and using the WCB metadata to select the WCB trajectories and extracting the"
	#"variables FOR THE ENTIRE SIMULATION. Returns them in simulation realtime and is NOT"
	#"reorganized in time. If a trajectory begins after the first rtime step, then the"
	#"variable is NaN until it starts."

	nan_offset = nest_offset_nan
	files = [glob("traj_sum*$(nf[end-15:end])",nest_path)[1] for nf in nest_WCB_fls]
	start_ind = nest_start 
	factor = nest_factor
	skp = nest_skp
	pth = nest_path

	#"Initialize array."
        v_tot = Vector{Matrix{Float64}}()
		
	#"Loop over all trajectory files."
        for f in files
		#"the first four time steps in the first couple of files are weird, thus"
		#"they are skipped."
		if f[end-18:end-8] == "tst00000003"
                        start_ind = 4
                elseif f[end-18:end-8] == "tst00001439"
			start_ind = 5
		else
                        start_ind = 1
                end
			
		#"Get data from trajectory file"
                ds = Dataset(f)

		#"Get metadata"
                w,a,t6 = WCB_info(f)

		#"Initialize array"
                v_f = fill(NaN,(length(w),192))

		#"Convert rtime to hours"
		r_t = ds["rtime"][start_ind:skp:end] ./ (3600*factor)

		#"Here variables are loaded in a loop that is multithreaded."
		#"If variable is NOT microphysical, take the ascent part and write to v_f array."
		#"If it is microphysical, it needs to go through the funciton "micro_singleline"."
		#"That function makes sure the variable is set to zero at beginning of ascent and"
		#"if fixes problems associated with ICON restarts. See the function for more info."
		if vr != "rtime"
                        vr_ar = ds[vr][:,start_ind:skp:end] .+ nan_offset
			@threads for i in 1:length(w)

			    frm1 = Int((r_t[2] - 0.5)*2)
			    frm2 = 1
			    
		 	    if vr ∉ microlist
                            	v_f[i,frm1:end] = vr_ar[w[i],:]
			    elseif vr in microlist
				v_f[i,frm1:end] = micro_singleline(r_t,vr_ar[w[i],:])
			    end	
			end
		elseif vr == "rtime"
			error("rtime isn't necessary in realtime")
		end
		v_tot = vcat(v_tot,[v_f[i,:] for i in 1:size(v_f)[1]])	
                close(ds)
        end
        return reduce(vcat,transpose.(v_tot))
end


function micro_singleline(rt_ar,vr_array)
	#"This function takes a 1-dim rtime and variable array and makes it such that"
	#"the variable is zero at the beginning of the array."
	#"This is for the microphysical variables which we want to start at zero at"
	#"the beginning of the ascent and accumulate with time."
	#"If the ICON simulation has to be restarted at a certain rtime, then the variable"
	#"will reset to zero in the simultion. To compensate for this, this function brings"
	#"the variable to the value it had before the restart. In the current form the function"
	#"can deal with 5 restarts. Just enter the rt values of the restart in ascending order."

	v = copy(vr_array)
	#"Restart times. -9999999 if no restart. rt5 MUST be smallest and rt1 largest."
	rt5 = 12.5
	rt4 = 48.5
	rt3 = 54.5
	rt2 = 66.5
	rt1 = -9999999
        
	ind1 = findall(rt_ar .== rt1)
        ind2 = findall(rt_ar .== rt2)
        ind3 = findall(rt_ar .== rt3)
        ind4 = findall(rt_ar .== rt4)
        ind5 = findall(rt_ar .== rt5)

        if length(ind1) .> 1
                error("ERROR1")
        elseif length(ind2) .> 1
                error("ERROR2")
        elseif length(ind3) .> 1
                error("ERROR3")
        elseif length(ind4) .> 1
                error("ERROR4")
        elseif length(ind5) .> 1
                error("ERROR5")
        end

	upto = length(v[(!isnan).(v)])
        if (length(ind1) > 0)
                if (ind1[1] != 1)
                        i = ind1[1]
                        v[i:upto] = v[i:upto] .+ v[i-1]
                end
        end
        if (length(ind2) > 0)
                if (ind2[1] != 1)
                        i = ind2[1]
                        v[i:upto] = v[i:upto] .+ v[i-1]
                end
        end
        if (length(ind3) > 0)
                if (ind3[1] != 1)
                        i = ind3[1]
                        v[i:upto] = v[i:upto] .+ v[i-1]
                end
        end
        if (length(ind4) > 0)
                if (ind4[1] != 1)
                        i = ind4[1]
                        v[i:upto] = v[i:upto] .+ v[i-1]
                end
        end
        if (length(ind5) > 0)
                if (ind5[1] != 1)
                        i = ind5[1]
                        v[i:upto] = v[i:upto] .+ v[i-1]
                end
        end
        return v
end

