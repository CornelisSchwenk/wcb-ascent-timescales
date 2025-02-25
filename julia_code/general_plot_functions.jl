#------------------- GENERAL PLOTTING FUNCTIONS -----------------
# "This script gets used by plot_functions.jl and it defines some general"
# "functions that can be used for plotting. icon_WCB_functions.jl must"
# "also be included so that these functions can work, because they use"
# "arrays or masks that are defined there."
#----------------------------------------------------------------

#--------"plot in the hours after ascent comparing conv and slow"---------
#"xar is the x-array with values being hours after the ascent."
#"vr_all is a dictionary with all variable values. vr_all[1] gives the"
#"values 1 hour after ascent and so on, and mxd is also a dictionary that"
#"gives the mask for which trajectories are still present x-hours after" 
#"the ascent. Both vr_all and mxd are generated in the respective plotting" 
#"function (i.e plot_Fig11 in plot_functions.jl)"

function plot_after_asc(xar,vr_all,mxd)
	p_out = plot(xar ./ 2,[nm.mean(vr_all[i][M_c[mxd[i]]]) for i in xar],label=false,
			 c=:blue,lw=2,xlabel="Hours after end of WCB ascent",
			 fillrange=([prcntl_1d(vr_all[i][M_c[mxd[i]]],10) for i in xar],
				    [prcntl_1d(vr_all[i][M_c[mxd[i]]],90) for i in xar]),
			 fillalpha=0.4)
	plot!(xar ./ 2,[nm.median(vr_all[i][M_c[mxd[i]]]) for i in xar],label=false,lw=1.5,
		  	c=:blue,ls=:dash)
	plot!(xar ./ 2,[nm.mean(vr_all[i][M_s[mxd[i]]]) for i in xar],label=false,c=:coral,lw=2,
		  fillrange=([prcntl_1d(vr_all[i][M_s[mxd[i]]],10) for i in xar],
				[prcntl_1d(vr_all[i][M_s[mxd[i]]],90) for i in xar]),fillalpha=0.4)
	plot!(xar ./ 2,[nm.median(vr_all[i][M_s[mxd[i]]]) for i in xar],label=false,lw=1.5,
		  c=:coral,ls=:dash)
	plot!([],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!([],c=:black,label="mean",lw=1.5)
        plot!([],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
        plot!([],c=:black,label="median",lw=1.5,ls=:dash)
        plot!(legend=:topright,legend_columns=2,legendfontsize=7)
	return p_out
end

#---"compare fast and convective trajs over normalised ascent time"--------
#"give this function an array that has the shape from arrays in dsnn"
#"and gives values along normalised ascent time and it will plot the"
#"evolution for slow and for convective trajectories with mean, median"
#"and 10th-90th percentiles shaded."
function plot_twotau_general(ar,av=false)
        if av==true
                a = ones(length(t6n),101)
                if !@isdefined rtav
                        rtav = dsnn["rtime"][:]
                        rtav = rtav .- rtav[:,1]
                end
                for i in 1:length(t6n)
                        a[i,2:end] = collect(xar_norm)[2:end] .* t6n[i]
                end
                ar = ar ./ a
        end
        fmc_10 = prcntl(ar[M_c,:],10)
        fmc_90 = prcntl(ar[M_c,:],90)

        fms_10 = prcntl(ar[M_s,:],10)
        fms_90 = prcntl(ar[M_s,:],90)

        p_out = plot(xar_norm,[nm.mean(ar[M_c,i]) for i in 1:101],
                   c=:blue,label="convective",lw=0,
                   xlabel="Normalised ascent time",
                   dpi=:300,fillrange=(fmc_10,fmc_90),fillalpha=0.5,
                   legendfontsize=10)

        plot!(xar_norm,[nm.median(ar[M_c,i]) for i in 1:101],c=:blue,label=false,lw=1,ls=:dash)
        plot!(xar_norm,[nm.mean(ar[M_c,i]) for i in 1:101],c=:blue,label=false,lw=1.2,ls=:solid)

        plot!(p_out,xar_norm,[nm.mean(ar[M_s,i]) for i in 1:101],lw=0,
                       c=:coral,label="slow",
                       fillrange=(fms_10,fms_90),fillalpha=0.5)

        plot!(xar_norm,[nm.median(ar[M_s,i]) for i in 1:101],c=:coral,label=false,lw=1,ls=:dash)
        plot!(xar_norm,[nm.mean(ar[M_s,i]) for i in 1:101],c=:coral,label=false,lw=1.2,ls=:solid)

        plot!([],lw=1.5,c=:black,label="mean")
        plot!([],lw=1.2,c=:black,ls=:dash,label="median")
        return p_out
end

#-----"normalize 2dim histogram"---
#"This function normalizes a 2dimension histogram along y."
function normalize_hist2d(hist,dim=1)
        zh = hist[1][1][:z]
        zh2 = replace(reshape(collect(zh),size(zh)),NaN => 0.0)
        norm = sum(zh2,dims=dim)
        normed = replace(zh2 ./ norm,0.0 => NaN)
        hist[1][1][:z] = Surface{Matrix{Float64}}(normed)
        hist[1][1][:clims_calculated] = (0.0,1.0)
end

#-----"general 2dim histogram over tau600"---
#"This function plots a 2d histogram of a 1dim array over tau_600."
function plot_hist2d_t6_general(ar,ybins,ylabel,loc=:auto,pl_mm=true)
	t6 = t6n
	h_ar = histogram2d((t6,ar),clims=(0,0.1),xlims=(0,48),ylims=(ybins[1],ybins[end]),
                         xlabel=taulabel,ylabel=ylabel,
                         bins=(0:1:48,ybins),labelfontsize=14,
                         colorbartitle="\n Normalised number of trajectories",
			 right_margin=5mm,left_margin=2mm,
			 legend=loc)
	normalize_hist2d(h_ar)
	if pl_mm == true
		plot!(1:0.5:48,t6_binning_fine(ar,t6,nm.mean),label="mean",lw=1.5,c=:cyan)
		plot!(1:0.5:48,t6_binning_fine(ar,t6,nm.median),label="median",
		      lw=1.5,c=:cyan,ls=:dash)
        end
	return h_ar
end


#-----"general 2dim histogram any two arrays"---
#"This function plots a 2d histogram of ar1 over ar2."
function plot_hist2d_2ars_general(ar1,ar2,xbins,ybins,xlabel,ylabel,loc=:auto,pl_mm=true)
	h_ar = histogram2d((ar1,ar2),clims=(0,0.1),xlims=(xbins[1],xbins[end]),
			   ylims=(ybins[1],ybins[end]),xlabel=xlabel,ylabel=ylabel,
			   bins=(xbins,ybins),labelfontsize=14,
			   colorbartitle="\n Normalised number of trajectories",
			   right_margin=5mm,legend=loc)
	normalize_hist2d(h_ar)
	if pl_mm == true
		plot!(xbins[2:end],bin_two_ars(ar2,ar1,xbins,nm.mean),
		       label="mean",lw=:2,c=:cyan)
		plot!(xbins[2:end],bin_two_ars(ar2,ar1,xbins,nm.median),
		       label="median",lw=:2,c=:cyan,ls=:dash)
	end
	return h_ar
end

#-------"Pressure bins"-------
#"This function bins an array vr into pressure bins. Vr must be from"
#"dsnr or dsnn. If masking is true it will take only values according"
#"to msk. Returns a dictionary with an array for each pressure bin, i.e"
#"P_d[500] has values for when p in (500,505]"
function get_var_Pbin(vr,taumasking=false,msk=false)
        y = vr
        if size(y)[2] == 192
                p = dsnr["p"][:] ./ 100
        elseif size(y)[2] == 101
                p = dsnn["p"][:] ./ 100
        else
                error("Take either from dsnr or from dsnn")
        end

        if taumasking==true
                y = y[msk,:]
                p = p[msk,:]
        end

        P_d = Dict()
        for pr in bin_ar
                inds = findall((p .> pr) .&& (p .<= pr+5))
                P_d[pr] = y[inds]
        end
        return P_d
end

#---------"Temperature bins"-----------
#"This function does the same but for temperature bins."
function get_var_Tbin(vr,taumasking=false,msk=false)
        y = vr

        if size(y)[2] == 192
                T = dsnr["t"][:] .- 273.15
        elseif size(y)[2] == 101
                T = dsnr["t"][:] .- 273.15
        else
                error("Take either from dsnr or from dsnn")
        end

        if taumasking==true
                y = y[msk,:]
                T = T[msk,:]
        end

        T_d = Dict()
        for temp in Tbin_ar
                inds = findall((T .> temp) .&& (T .<= temp+1))
                T_d[temp] = y[inds]
        end
        return T_d
end

#----------"Bin arithmetics"-------
#"This function returns the prc-percentile of a binned dictionary"
#"obtained from either get_var_Pbin or get_var_Tbin."
function bin_prcntl(bindict,prc)
        #"This function takes a dictionary created by get_var_pbin (that function bins WCB data"
        #"from the WCB_vars... files into pressure bins). Then it gives you the prc percentile."
        #"The keys are the pressure bins."
        x_ar = sort([k for k in keys(bindict)])
        prc_ar = []
        for x in x_ar
                l = length(bindict[x])/100
                ind = Int(round(l*prc,digits=0))
                push!(prc_ar,sort(bindict[x])[ind])
        end
        return prc_ar
end

#"This function returns an array with the mean for each P or T bin."
function bin_mean(bindict)
        #"Takes the mean for each pressure bin and can accomodate NaNs"
        x_ar = sort([k for k in keys(bindict)])
        return [nm.mean(bindict[x]) for x in x_ar]
end

#"This function returns an array with the median for each P or T bin."
function bin_median(bindict)
        #"Takes the median for each pressure bin and can accomodate NaNs"
        x_ar = sort([k for k in keys(bindict)])
        return [nm.median(bindict[x]) for x in x_ar]
end

#"Plot and array in temperature bins. vr must be from dsnn or dsnr or"
#"at least derived from them (i.e vr = dsnn["qc"][:] .+ dsnn["qr"][:])"
function plot_Tbinned_var_twotau(vr,fix=false)
        vn1 = get_var_Tbin(vr,true,M_c)
        vn2 = get_var_Tbin(vr,true,M_s)

        vn1_median = bin_median(vn1)
        vn1_mean = bin_mean(vn1)
        f1_10 = bin_prcntl(vn1,10) .- 1e-20
        f1_90 = bin_prcntl(vn1,90) .+ 1e-20

        vn2_median = bin_median(vn2)
        vn2_mean = bin_mean(vn2)
        f2_10 = bin_prcntl(vn2,10) .- 1e-20
        f2_90 = bin_prcntl(vn2,90) .+ 1e-20

        if fix == true
                f1_10[f1_10 .<= 0] .= 0.1
                f2_10[f2_10 .<= 0] .= 0.1
        end

        pln = plot(Tbin_ar,vn1_mean,xflip=false,lw=0,c=:blue,
                   fillrange=(f1_10,f1_90),
                   fillalpha=0.4,xlabel=L"\mathrm{Temperature}\quad[{}^{\circ}\mathrm{C}]",
                   label="convective",legend=:topright,
                   labelfontsize=14,legendfontsize=10,dpi=:400)

        plot!(pln,Tbin_ar,vn2_mean,xflip=false,lw=0,c=:coral,
                   fillrange=(f2_10,f2_90),
                   fillalpha=0.4,label="slow")

        plot!(pln,Tbin_ar,vn1_mean,lw=1.5,ls=:solid,c=:blue,label=false)
        plot!(pln,Tbin_ar,vn2_mean,lw=1.5,ls=:solid,c=:coral,label=false)

        plot!(pln,Tbin_ar,vn1_median,lw=1.5,ls=:dash,c=:blue,label=false)
        plot!(pln,Tbin_ar,vn2_median,lw=1.5,ls=:dash,c=:coral,label=false)

        plot!([],lw=1.5,c=:black,label="mean")
        plot!([],lw=1.2,c=:black,ls=:dash,label="median")
	return pln
end

#"Plot and array in pressure bins. vr must be from dsnn or dsnr or at"
#"least derived from them (i.e vr = dsnn["qc"][:] .+ dsnn["qr"][:])"
function plot_Pbinned_var_twotau(vr)
        vn1 = get_var_Pbin(vr,true,M_c)
        vn2 = get_var_Pbin(vr,true,M_s)

        vn1_median = bin_median(vn1)
        vn1_mean = bin_mean(vn1)
        f1_10 = bin_prcntl(vn1,10) .- 1e-20
        f1_90 = bin_prcntl(vn1,90) .+ 1e-20

        vn2_median = bin_median(vn2)
        vn2_mean = bin_mean(vn2)
        f2_10 = bin_prcntl(vn2,10) .- 1e-20
        f2_90 = bin_prcntl(vn2,90) .+ 1e-20

        pln = plot(bin_ar,vn1_mean,xflip=true,lw=0,c=:blue,
                   fillrange=(f1_10,f1_90),
                   fillalpha=0.4,xlabel=L"\mathrm{Pressure}\quad[\mathrm{hPa}]",
                   label="convective",legend=:topright,
                   labelfontsize=14,legendfontsize=10,dpi=:400)

        plot!(pln,bin_ar,vn2_mean,xflip=true,lw=0,c=:coral,
                   fillrange=(f2_10,f2_90),
                   fillalpha=0.4,label="slow")
        
	plot!(pln,bin_ar,vn1_mean,lw=1.5,ls=:solid,c=:blue,label=false)
        plot!(pln,bin_ar,vn2_mean,lw=1.5,ls=:solid,c=:coral,label=false)

        plot!(pln,bin_ar,vn1_median,lw=1.5,ls=:dash,c=:blue,label=false)
        plot!(pln,bin_ar,vn2_median,lw=1.5,ls=:dash,c=:coral,label=false)
	
	plot!([],lw=1.5,c=:black,label="mean")
	plot!([],lw=1.2,c=:black,ls=:dash,label="median")
	return pln
end
