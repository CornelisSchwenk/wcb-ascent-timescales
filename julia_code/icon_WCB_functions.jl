#--------------------------------------------------------------------------
#       "Skript for general work with ICON WCB trajectory files etc"
#       "Created by Cornelis Schwenk 27.11.2023"
#       "c.schwenk@uni-mainz.de"

#	"This script should only be included once all of the WCB_vars... files"
#	"have been created. See WCB_selection.jl"
#---------------------------------------------------------------------------

include("use_config.jl")
using Glob, NCDatasets, Plots
using Plots.PlotMeasures
using Statistics
using BenchmarkTools, Dates, DataStructures
using DelimitedFiles
using NaNMath; nm=NaNMath
using StatsBase
using LaTeXStrings
using ColorSchemeTools
using LsqFit

include("thermodynamic_functions.jl")
include("statistical_functions.jl")
include("WCB_selection.jl")

#"First we define some files and arrays we can use generally."
#"Paths and so on are defined in WCB_selection.jl"

#----"files from Eulerian ICON output (empty if not available)"------
nest_NWP_fls_dom1 = glob(nest_NWP_path*"*mod_NWP*DOM01*.nc")
nest_NWP_fls_dom2 = glob(nest_NWP_path*"*mod_NWP*DOM02*.nc")
nest_NWP_fls_dom3 = glob(nest_NWP_path*"*mod_NWP*DOM03*.nc")

#----"define arrays from Eulerian data if available"----
if (length(nest_NWP_fls_dom1) != 0) .& (length(nest_NWP_fls_dom2) != 0) .& (length(nest_NWP_fls_dom3) != 0)
    nest_dom1_ds_NWP = Dataset(nest_NWP_fls_dom1[1])
    nest_NWP_lon_dom1 = Array(nest_dom1_ds_NWP["lon"])
    nest_NWP_lat_dom1 = Array(nest_dom1_ds_NWP["lat"])

    nest_dom2_ds_NWP = Dataset(nest_NWP_fls_dom2[1])
    nest_NWP_lon_dom2 = Array(nest_dom2_ds_NWP["lon"])
    nest_NWP_lat_dom2 = Array(nest_dom2_ds_NWP["lat"])

    nest_dom3_ds_NWP = Dataset(nest_NWP_fls_dom3[1])
    nest_NWP_lon_dom3 = Array(nest_dom3_ds_NWP["lon"])
    nest_NWP_lat_dom3 = Array(nest_dom3_ds_NWP["lat"])

    #	"Domain borders"
    global maxlon1 = maximum(nest_NWP_lon_dom1)
    global minlon1 = minimum(nest_NWP_lon_dom1)
    global maxlat1 = maximum(nest_NWP_lat_dom1)
    global minlat1 = minimum(nest_NWP_lat_dom1)

    global maxlon2 = maximum(nest_NWP_lon_dom2)
    global minlon2 = minimum(nest_NWP_lon_dom2)
    global maxlat2 = maximum(nest_NWP_lat_dom2)
    global minlat2 = minimum(nest_NWP_lat_dom2)

    global maxlon3 = maximum(nest_NWP_lon_dom3)
    global minlon3 = minimum(nest_NWP_lon_dom3)
    global maxlat3 = maximum(nest_NWP_lat_dom3)
    global minlat3 = minimum(nest_NWP_lat_dom3)
end

#----"Define arrays used for plotting"----
global bin_ar = collect(240:5:845)
global Tbin_ar = collect(-40:1:20)
global xar_norm = LinRange(0,1,101)

#----"Load specific arrays if the proper files are available"----
if length(glob(nest_path*"WCB_tau_vars_realtime.nc")) > 0
    global ds_nest_WCB_rt = Dataset(nest_path*"WCB_tau_vars_realtime.nc")
    global dsnr = ds_nest_WCB_rt
    if (@isdefined pnr) == false
        global pnr = Array(dsnr["p"]) ./ 100
    end
    if !@isdefined t6n
        global t6n = Array(dsnr["tau_600"])
    end
    if !@isdefined dp2
        global dp2 = Array(dsnr["dp2"])
    end
else
    error("If you don't have the WCB_tau_vars_realtime.nc file, please
        download it from zenodo and place it in ../nest_data")
end

if length(glob(nest_path*"WCB_tau_vars_normed.nc")) > 0
    global ds_nest_WCB_normed = Dataset(nest_path*"WCB_tau_vars_normed.nc")
    global dsnn = ds_nest_WCB_normed
else
        println("----- Ignore this message if you are currently executing the 
        function write_tauWCB_file_normed() ------ 
        Otherwise: you seem to NOT have the WCB_tau_vars_normed.nc file, 
        so you have to first download WCB_tau_vars_realtime.nc from zenodo, 
        place it in ../nest_data, and then create the normed file using 
        the function write_tauWCB_var_file_normed() from WCB_selection.jl")
end

if ((@isdefined wn) == false) .& ((@isdefined pnr) == true) .& (length(nest_fls) > 0)
    #"wn: wcb indices (in y) in traj files."
    #"an: indices (in x) of time when ascent begins."
    #"tn: tau_600 in time steps."
    #"ns: indices (in x) the NaNs end in the _realtime file."
    #"ains: CartesianIndices of all ascent starting points"
    global wn,an,tn = get_w_a_all()
    global ns = [findfirst((!isnan).(pnr[i,:]))-1 for i in 1:length(tn)]
    global ains = [CartesianIndex(i,an[i] + ns[i]) for i in 1:length(tn)]
end

GC.gc()


#---- Cloud parameters-----
global a_geo_ice = 0.835
global b_geo_ice = 0.390

#"Define global masks for convective trajectories (M_c) and"
#"slow trajectories (M_s)"

global M_c = t6n .< 5
global M_s = (t6n .> 20) .& (dp2 .< 350)
global M_n = 18 .>= t6n .>= 6


#------------------------------------------------------------
#		Functions
#------------------------------------------------------------

#	"Conversion of rtime to DateTime and vice versa"
rt_to_dt(rt) = DateTime(2017,9,20,0) + Minute(rt*60)
dt_to_rt(dt::DateTime) = Minute(dt - DateTime(2017,9,20,0)).value / 60

function find_max_dp2()
    #"This function finds the fastest 2hr pressure ascent"
    #"during the WCB-ascent"
    p_ar = pnr
    if !(@isdefined ains)
        error("ains needs to be defined. For this you need to have downloaded the ICON data, especially the trajectories, and run through the WCB selection process to create the merged and summated files as well as the WCB_tau_tst* files. If you need dp2, take it from the WCB_tau_vars_realtime.nc file that is provided in zenodo!")
    end
    ain = ains	
    
    dp2_out = zeros(1:size(p_ar)[1])
    for k in 1:length(dp2_out)
        p_k = p_ar[k,ain[k][2]:end]
        dp2_1 = maximum(p_k[1:end-4] .- p_k[5:end])
        dp2_2 = maximum(p_k[2:end-4] .- p_k[6:end])
        dp2_3 = maximum(p_k[3:end-4] .- p_k[7:end])
        dp2_4 = maximum(p_k[4:end-4] .- p_k[8:end])
        dp2_out[k] = maximum([dp2_1,dp2_2,dp2_3,dp2_4])
    end
    return dp2_out
end

function find_CAP()
    #"This function finds the fastest 2hr pressure ascent"
    #"during the WCB-ascent"
    p_ar = pnr
    if !(@isdefined ains)
        error("ains needs to be defined. For this you need to have downloaded the ICON data, especially the trajectories, and run through the WCB selection process to create the merged and summated files as well as the WCB_tau_tst* files. If you need dp2, take it from the WCB_tau_vars_realtime.nc file that is provided in zenodo!")
    end
    ain = ains	
    
    CAP_out = zeros(1:size(p_ar)[1])
    for k in 1:length(CAP_out)
        p_k = p_ar[k,ain[k][2]:end]
        CAP_out[k] = find_max_dpx(p_k,1)
    end
    return CAP_out
end

function calc_D_ice(q,qn)
    m = q / qn
    return a_geo_ice * (m ^ b_geo_ice)
end

function calc_radius_ice(q,qn)
    m = q ./ qn
    return 0.5 .* a_geo_ice .* (m .^ b_geo_ice) 
end

function calc_fx_ice(q,qn)
    x = q ./ qn
    return 10 .* qn .* exp( 60 .^ (1/3) ) ./ x 
end

function calc_vertical_velocity()
    #"returns ascent velocity in m/s derived from altitude change over runtime"
    ds = dsnn
    alt = Array(ds["alt"])
    rt = Array(ds["rtime"])
    altd = alt[:,2:end] .- alt[:,1:end-1]
    rtd = rt[:,2:end] .- rt[:,1:end-1]
    wv = altd ./ (rtd .* 60 .* 60)
    return wv
end

function calc_vertical_velocity_realtime()
    #"returns ascent velocity in m/s derived from altitude change over runtime"
    ds = dsnr
    alt = Array(ds["alt"])
    altd = alt[:,2:end] .- alt[:,1:end-1]
    rtd = 0.5
    wv = altd ./ (rtd .* 60 .* 60)
    return wv
end

function t6_binning(ar,t6,func)
    #"Bin func(values) into tau_600 bins."
    #"Used to calculate median and mean as function of t6."
    #"Example: CR_mean = t6_binning(CR,tau_600,mean)."
    #"Uses bin_two_ars function from statistical_functions.jl"
    bins = 1:48
    return bin_two_ars(ar,t6,bins,func)
end

function t6_binning_fine(ar,t6,func)
        v = []
        for t in collect(1:0.5:48)
                if length(ar[t-1 .< t6 .<= t]) == 0
                        push!(v,NaN)
                else
                        push!(v,func(ar[t-1 .< t6 .<= t]))
                end
        end
        return v
end

function t6_prcntl(ar,t6,prc)
    v = Float64[]
    for t in collect(1:0.5:48)
        if length(ar[t-1 .< t6 .<= t]) == 0
            push!(v,NaN)
        else
            s_ar = sort(ar[t-1 .< t6 .<= t])
            iprc = Int(round(length(s_ar)*prc/100,digits=0))
            push!(v,s_ar[iprc])
        end
    end
    return v
end

function anyar_prcntl(ar1,ar2,bins,prc)
    #ar1 is binned into bins of ar2
    v = Float64[]
    for i in 1:(length(bins)-1)
        if length(ar1[bins[i] .< ar2 .<= bins[i+1]]) < 100
            push!(v,NaN)
        else
            s_ar = sort(ar1[bins[i] .< ar2 .<= bins[i+1]])
            push!(v,prcntl_1d(s_ar,prc))
        end
    end
    return v
end


function arnorm(ar)
    return (ar .- mean(ar)) ./ std(ar)
end

function taumask(a1,a2)
    #"Create a mask for a1 .< tau_600 .< a2"
        tau = Array(dsnn["tau_600"]) ./ 2
        return (tau .> a1) .&& (tau .<= a2)
end


function find_min_dpx(p,x)
    #"takes 1-dim p_array and time steps x"
    #"finds the minimum ascent in x time steps"
    dpx = p[1:end-x] .- p[x+1:end]
    return minimum(dpx)
end

function find_max_dpx(p,x)
    #"takes 1-dim p_array and time steps x"
    #"finds the maximum ascent in x time steps"
    dpx = p[1:end-x] .- p[x+1:end]
    return maximum(dpx)
end

function flux(hy::String,form::String)
    #"Calculate precipitation flux of specific species"
    #"Chose normed, regular or realtime to specify which WCB_tau_vars.... to use"
    if form == "normed"
        ds = dsnn
    elseif form =="realtime"
        ds = dsnr
    else
        error("form must be 'realtime' or 'normed'")
    end
    in = Array(ds["q$(hy)_in"])
    out = Array(ds["q$(hy)_out"])
    
    return in .- out
end

function flux_fast_ind(hy::String,form::String,ind::Int)
    #"Faster version for specific time step"
    if form == "normed"
        ds = dsnn
    elseif form =="realtime"
        ds = dsnr
    else
        error("form must be 'realtime' or 'normed'")
    end
    in = ds["q$(hy)_in"][:,ind]
    out = ds["q$(hy)_out"][:,ind]

    return in .- out
end

function flux_superquick_ind(hy::String,ds,ind::Int)
    #"Faster version for specific time step, input ds"
    in = ds["q$(hy)_in"][:,ind]
    out = ds["q$(hy)_out"][:,ind]
    return in .- out
end

#	"Faster version for just the last time step"
function flux_fast_last(hy::String,form::String)
    if form == "normed"
        ds = dsnn
    elseif form =="realtime"
        ds = dsnr
    else
        error("form must be 'realtime' or 'normed'")
    end
    
    in = ds["q$(hy)_in"][:,end]
    out = ds["q$(hy)_out"][:,end]

    return in .- out
end

function mean_flux(hy::String,form::String)
    #"Calculate the mean across time axis"
    f = flux(hy,form)	
    return [nm.mean(f[:,i]) for i in 1:size(f)[2]]
end

function get_tWCB_plus_ts_indices(ts::Int)
    #"This function gives the CartesianIndices of trajectories in the"
    #"realtime file ts-timestepts after the end of WCB ascent. The msk" 
    #"is a mask that is false if ts-hours is outside of the sim."
    ds = dsnr
    upto = length(Array(dsnr["tau_600"]))
    ts_inds = [CartesianIndex(k,find_stop(k) + ts) for k in 1:upto]
    msk = [ts_inds[k][2] for k in 1:upto] .<= 192
    return (ts_inds,msk)
end

#------------------------------------------------
#	"Calculations for water budget and residuals etc."
#	"See publication for details."
#	"The following functions only use the time-normalized data from WCB_tau_vars_normed.nc"
#------------------------------------------------

function calc_P()
        return sum([flux(hy,"normed") for hy in hymet[2:end]])
end

function calc_dqv_traj(ind=101::Int)
        #"calculate change in moisture along ascent"
        ds = dsnn
		dqv_traj = -ds["qihh"][:,ind] .- ds["qc_nuc"][:,ind] .- ds["cond"][:,ind] .-
					ds["evap"][:,ind] .+ ds["r_evap"][:,ind] .+ ds["f_evap"][:,ind] .-
					ds["qx_dep"][:,ind] .- ds["satad2"][:,ind] .+ ds["qvturb"][:,ind] .+
					ds["qvconv"][:,ind]
        return dqv_traj
end

function calc_dqv_traj_all(ds)
        #"calculate change in moisture along ascent"
		dqv_traj = -Array(ds["qihh"]) .- Array(ds["qc_nuc"]) .- Array(ds["cond"]) .-
					Array(ds["evap"]) .+ Array(ds["r_evap"]) .+ Array(ds["f_evap"]) .-
					Array(ds["qx_dep"]) .- Array(ds["satad2"]) .+ Array(ds["qvturb"]) .+
					Array(ds["qvconv"])
        return dqv_traj
end

function calc_dqhy_traj(ind=101::Int)
        #"calculate change in hydrometeors along ascent"
        ds = dsnn
        allfl = sum([flux_fast_ind(hy,"normed",ind) for hy in hymet[2:end]])
		dqhy_traj = ds["qihh"][:,ind] .+ ds["qc_nuc"][:,ind] .+ ds["cond"][:,ind] .+
				ds["evap"][:,ind] .-
				ds["r_evap"][:,ind] .- ds["f_evap"][:,ind] .+ ds["qx_dep"][:,ind] .+
				ds["satad2"][:,ind] .+ ds["qcturb"][:,ind] .+ ds["qcconv"][:,ind] .+
				ds["qsconv"][:,ind] .+ ds["qiconv"][:,ind] .+ ds["qrconv"][:,ind] .- allfl
		return dqhy_traj
	end

	function calc_dqhy_traj_superquick(ds,ind::Int)
		#"calculate change in hydrometeors along ascent"
		allfl = sum([flux_superquick_ind(hy,ds,ind) for hy in hymet[2:end]])
		dqhy_traj = ds["qihh"][:,ind] .+ ds["qc_nuc"][:,ind] .+ ds["cond"][:,ind] .+
				ds["evap"][:,ind] .-
				ds["r_evap"][:,ind] .- ds["f_evap"][:,ind] .+ ds["qx_dep"][:,ind] .+
				ds["satad2"][:,ind] .+ ds["qcturb"][:,ind] .+ ds["qcconv"][:,ind] .+
				ds["qsconv"][:,ind] .+ ds["qiconv"][:,ind] .+ ds["qrconv"][:,ind] .- allfl
        return dqhy_traj
end

function calc_dqhy_traj_all(ds,allfl)
        #"calculate change in hydrometeors along ascent"
	dqhy_traj = Array(ds["qihh"]) .+ Array(ds["qc_nuc"]) .+ Array(ds["cond"]) .+ Array(ds["evap"]) .-
				Array(ds["r_evap"]) .- Array(ds["f_evap"]) .+ Array(ds["qx_dep"]) .+
				Array(ds["satad2"]) .+ Array(ds["qcturb"]) .+ Array(ds["qcconv"]) .+
				Array(ds["qsconv"]) .+ Array(ds["qiconv"]) .+ Array(ds["qrconv"]) .- allfl
	return dqhy_traj
end

function calc_res_qv(ind=101::Int)
        #"calculate residual water vapor along ascent"
        ds = dsnn
        res_qv = ds["qv"][:,1] .- ds["qv"][:,ind] .+ calc_dqv_traj(ind)
        return res_qv
end

function calc_res_qv_all(ds)
        #"calculate residual water vapor along ascent"
        res_qv = ds["qv"][:,1] .- Array(ds["qv"]) .+ calc_dqv_traj_all(ds)
        return res_qv
end

function calc_res_qhy(ind=101::Int)
        #"calculate residual hydrometeors along ascent"
        ds = dsnn
        hy_tot_ind = sum([ds["q$(hy)"][:,ind] for hy in hymet])
        qhy0 = sum([ds["q$(hy)"][:,1] for hy in hymet])
        res_qhy = qhy0 .- hy_tot_ind .+ calc_dqhy_traj(ind)
        return res_qhy
end

function calc_res_qhy_quick(qhy0,ind=101::Int)
        #"calculate residual hydrometeors along ascent"
        ds = dsnn
        hy_tot_ind = sum([ds["q$(hy)"][:,ind] for hy in hymet])
        res_qhy = qhy0 .- hy_tot_ind .+ calc_dqhy_traj(ind)
        return res_qhy
end

function calc_res_qhy_superquick(ds,qhy0,ind::Int)
        #"calculate residual hydrometeors along ascent"
        hy_tot_ind = sum([ds["q$(hy)"][:,ind] for hy in hymet])
        res_qhy = qhy0 .- hy_tot_ind .+ calc_dqhy_traj_superquick(ds,ind)
        return res_qhy
end

function calc_res_qhy_all(ds,qhy,allfl)
        #"calculate residual hydrometeors along ascent"
        res_qhy = qhy[:,1] .- qhy .+ calc_dqhy_traj_all(ds,allfl)
        return res_qhy
end

function calc_res(ind=101::Int)
        return calc_res_qhy(ind) .+ calc_res_qv(ind)
end

function calc_res_all(ds,qhy,allfl)
        return calc_res_qhy_all(ds,qhy,allfl) .+ calc_res_qv_all(ds)
end

function calc_HYD(ind=101::Int)
        ds = dsnn
		return sum([ds["q$(hy)"][:,1] for hy in hymet]) .+ ds["qcturc"][:,ind] .+
				ds["qcconv"][:,ind] .+ ds["qrconv"][:,ind] .+ ds["qsconv"][:,ind] .+
				ds["qiconv"][:,ind] .- calc_res_qhy(ind)
	end

function calc_HYD_quick(qhy0,ind=101)
	ds = dsnn
	return qhy0 .+ ds["qcturc"][:,ind] .+ ds["qcconv"][:,ind] .+ ds["qrconv"][:,ind] .+ 
		ds["qsconv"][:,ind] .+ ds["qiconv"][:,ind] .- 
		calc_res_qhy_quick(qhy0,ind)
end

function calc_HYD_superquick(ds,qhy0,ind::Int)
	return qhy0 .+ ds["qcturc"][:,ind] .+ ds["qcconv"][:,ind] .+ ds["qrconv"][:,ind] .+ 
		ds["qsconv"][:,ind] .+ ds["qiconv"][:,ind] .- 
		calc_res_qhy_superquick(ds,qhy0,ind)
end
end

function calc_HYD_all(ds,qhy,allfl)
	one = qhy[:,1] .+ Array(ds["qcturc"]) .+ Array(ds["qcconv"]) .+ Array(ds["qrconv"])
	two = Array(ds["qsconv"]) .+ Array(ds["qiconv"])
	GC.gc()
	return one .+ two .- calc_res_qhy_all(ds,qhy,allfl)
end

function calc_VAP(ind=101::Int)
	ds = dsnn
	return ds["qv"][:,1] .+ ds["qvturc"][:,ind] .+ ds["qvconv"][:,ind] .- calc_res_qv(ind)
end

function calc_VAP_all(ds)
	return ds["qv"][:,1] .+ Array(ds["qvturc"]) .+ Array(ds["qvconv"]) .- calc_res_qv_all(ds)
end

function calc_qtot(ind=101::Int)
	ds = dsnn
	return sum([ds["q$(hy)"][:,ind] for hy in ["v","c","r","i","s","g","h"]])
end

function calc_qtot_all()
        ds = dsnn
	qtot = Array(ds["qv"])
	for hy in ["c","r","i","s","g","h"]
		qtot = qtot .+ Array(ds["q$(hy)"])
	end
	GC.gc()
	return qtot
end

function calc_convs(ind=101::Int)
        ds = dsnn
        return sum([ds["q$(x)conv"][:,ind] for x in ["v","c","r","i","s"]])
end

function calc_convs_all(ds)
        convs = Array(ds["qvconv"])
	for hy in ["c","r","i","s"]
		convs = convs .+ Array(ds["q$(hy)conv"])
	end
	return convs
end

function calc_turbs(ind=101::Int)
        ds = dsnn
        return sum([ds["q$(hy)turb"][:,ind] for hy in ["c","v"]])
end

function calc_turbs_all(ds)
        return sum([Array(ds["q$(hy)turb"]) for hy in ["c","v"]])
end

function calc_turcs(ind=101::Int)
        ds = dsnn
        return sum([ds["q$(hy)turc"][:,ind] for hy in ["c","v"]])
end

function calc_qtc(ind=101::Int)
        return calc_convs(ind) .+ calc_turbs(ind)
end

function calc_qtc_all(ds)
        return calc_convs_all(ds) .+ calc_turbs_all(ds)
end

function calc_qtcr(ind=101::Int)
        return calc_res(ind) .- calc_qtc(ind)
end

function calc_qtcr_all(ds,qhy,allfl)
        return calc_res_all(ds,qhy,allfl) .- calc_qtc_all(ds)
end

function calc_C_hy(ind=101::Int)
	ds = dsnn
	C_hy = ds["cond"][:,ind] .+ ds["qc_nuc"][:,ind] .+ ds["qihh"][:,ind] .+ ds["depo"][:,ind]
	return C_hy
end

function calc_E_v(ind=101::Int)
	ds = dsnn
	E_v = ds["evap"][:,ind] .- ds["r_evap"][:,ind] .- ds["f_evap"][:,ind] .+ ds["subl"][:,ind]
	return E_v
end

function calc_eps(ind=101::Int)
        return calc_HYD(ind) ./ calc_VAP(ind)
end

function calc_PE(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
	ds = dsnn
        tot_prec = sum([flux_fast_ind(i,"normed",ind) for i in hymet[2:end]])
        qhy0 = sum([ds["q$(hy)"][:,1] for hy in hymet])
	HYD = calc_HYD_quick(qhy0,ind)
        PE = tot_prec ./ (C_hy .+ E_v .+ HYD)
	PE[PE .> 1] .= NaN
	PE[PE .< 0] .= NaN
	return PE
end

function calc_PE_time(qhy,allfl)
	ds = dsnn
	
	C_hy = Array(ds["cond"]) .+ Array(ds["qc_nuc"]) .+ Array(ds["qihh"]) .+ Array(ds["depo"])
	E_v = Array(ds["evap"]) .- Array(ds["r_evap"]) .- Array(ds["f_evap"]) .+ Array(ds["subl"]) 
		
	GC.gc()
	HYD = calc_HYD_all(ds,qhy,allfl)

	PE = allfl ./ (C_hy .+ E_v .+ HYD)
	PE[PE .> 1] .= NaN
	PE[PE .< 0] .= NaN

	GC.gc()
	return PE
end

function calc_PE_Eul(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        P1 = sum([flux_fast_ind(i,"normed",ind) for i in hymet[2:end]])
        PE = P1 ./ (C_hy .+ E_v)
        return PE
end

function calc_PE_notc(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        tot_prec = sum([flux_fast_ind(i,"normed",ind) for i in hymet[2:end]])
        qhy0 = sum([dsnn["q$(hy)"][:,1] for hy in hymet])
        return tot_prec ./ (C_hy .+ E_v .+ qhy0)
end

function calc_CR(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        VAP = calc_VAP(ind)
        CR = (C_hy .+ E_v) ./ VAP
	CR[CR .> 1] .= NaN
	CR[CR .< 0] .= NaN
	return CR
end

function calc_CR_time()
	ds = dsnn
		C_hy = Array(ds["cond"]) .+ Array(ds["qc_nuc"]) .+ Array(ds["qihh"]) .+ Array(ds["depo"])
		E_v = Array(ds["evap"]) .- Array(ds["r_evap"]) .- Array(ds["f_evap"]) .+ Array(ds["subl"]) 
		VAP = calc_VAP_all(ds)
		CR = (C_hy .+ E_v) ./ VAP
	CR[CR .> 1] .= NaN
	CR[CR .< 0] .= NaN
	return CR
end

function calc_CR_Eul(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        Qtot0 = calc_qtot(1)
        CR = (C_hy .+ E_v) ./ Qtot0
        return CR
end

function calc_CR_notc(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        qv0 = dsnn["qv"][:,1]
        return (C_hy .+ E_v) ./ qv0
end


function calc_DR(ind=101::Int)
        qtot_ind = calc_qtot(ind)
        qtot_1 = calc_qtot(1)
        return (qtot_1 .- qtot_ind) ./ (qtot_1)
end

function calc_DR_time()
        qtot_t = mapreduce(permutedims,vcat,[calc_qtot(i) for i in 1:101])'
        qtot_1 = calc_qtot(1)
        return (qtot_1 .- qtot_t) ./ (qtot_1)
end

function calc_DR_mix(ind=101::Int)
        qtot0 = calc_qtot(1)
        qtcr = calc_qtcr(ind)
        return qtcr ./ qtot0
end

function calc_DR_mix_time(qhy,allfl)
	ds = dsnn
        qtot0 = calc_qtot(1)
	qtcr = calc_qtcr_all(ds,qhy,allfl)
        return qtcr ./ qtot0
end

function calc_DR_mphys(ind=101::Int)
        P1 = sum([flux_fast_ind(hy,"normed",ind) for hy in hymet[2:end]])
        qtot0 = calc_qtot(1)
        qtcr = calc_qtcr(ind)
	DR_mphys =  P1 ./ (qtot0 .- qtcr)
	DR_mphys[DR_mphys .> 1] .= NaN
	DR_mphys[DR_mphys .< 0] .= NaN
	return DR_mphys
end

function calc_DR_mphys_time(qhy,allfl)
	ds = dsnn
	qtot0 = calc_qtot(1)
	qtcr = calc_qtcr_all(ds,qhy,allfl)
	DR_mphys =  allfl ./ (qtot0 .- qtcr)
	DR_mphys[DR_mphys .> 1] .= NaN
	DR_mphys[DR_mphys .< 0] .= NaN
	return DR_mphys
end

function calc_DR_res(ind=101::Int)
        resqv = calc_res_qv(ind)
        resqhy = calc_res_qhy(ind)
        qtot = calc_qtot(1)
        return (resqv .+ resqhy) ./ qtot

end

function calc_q_removal(ind=101::Int)
        #"this function tracks the changes in qtot budget"
        #"the last time-entry is equal to change in qtot"
        ds = dsnn
        allfl = sum([flux_fast_ind(hy,"normed",ind) for hy in hymet[2:end]])
        turbs = sum([ds["q$(hy)turb"][:,ind] for hy in ["c","v"]])
        convs = sum([ds["q$(x)conv"][:,ind] for x in ["v","c","r","i","s"]])
        resqv = calc_res_qv(ind)
        resqhy = calc_res_qhy(ind)
        return allfl .- turbs .- convs .+ resqv .+ resqhy
end

function calc_delta_H()
        ds = dsnn
        return Array(ds["qihh"]) .+ Array(ds["qc_nuc"]) .+ Array(ds["cond"]) .+ Array(ds["evap"]) .+ Array(ds["r_evap"]) .+
                Array(ds["f_evap"]) .+ Array(ds["qx_dep"]) .+ Array(ds["satad2"])
end
