#----------------------------------------------------------------------------
#	"This code defines the functions to make most of the plots for the paper"
#	"Cornelis Schwenk 09.07.2024"
#	"c.schwenk@uni-mainz.de"
#----------------------------------------------------------------------------

include("icon_WCB_functions.jl")
using StatsPlots
include("general_plot_functions.jl")

savepath = "../plots/"

#"Math label used for plotting" 
taulabel = L"\tau_{600}\quad [\mathrm{h}]" 

if length(glob("../nest_data/PECR.nc")) == 0
	error("You need to download the file PECR.nc from zenodo or create it using write_PECR.jl")
end

dsPE = Dataset("../nest_data/PECR.nc");

#"Figure 4: tau600 and tau400 hist"
function plot_Fig4(save=false)
	p1 = plot(t6n,bins=(0:0.5:48),st=:stephist,
                      label=L"\tau_{600}",xlabel=L"\mathrm{Time}\quad[\mathrm{h}]",
                      ylabel=L"\mathrm{Counts}",dpi=:500,
                      lw=2,c=:black,labelfontsize=14,legendfontsize=10,
		      title="(a)",titlefontsize=10,title_location=:left)
        tWCB = Array(dsnn["tau_WCB"])
        plot!(t6n[M_c],bins=(0:0.5:48),st=:stephist,
               fill=:true,c=:blue,alpha=:0.5,lw=0,label=L"\mathrm{convective}")
        plot!(t6n[M_n],bins=(0:0.5:48),st=:stephist,
               fill=:true,c=:grey,alpha=:0.5,lw=0,label=L"\mathrm{normal}")
        plot!(t6n[M_s],bins=(0:0.5:48),st=:stephist,
               fill=:true,c=:coral,alpha=:0.5,lw=0,label=L"\mathrm{slow}")
        plot!(tWCB,bins=(0:0.5:55),st=:stephist,label=L"\tau_\mathrm{WCB}",lw=1.5,
               c=:black,ls=:dot)

	t3n = [find_tau_x00(pnr[k,an[k]+ns[k]:end],300) for k in 1:length(t6n)] ./ 2
	t4n = [find_tau_x00(pnr[k,an[k]+ns[k]:end],400) for k in 1:length(t6n)] ./ 2
	p2 = plot(t3n,bins=(0:0.5:30),st=:stephist,
                      label=L"\tau_{300}",xlabel=L"\mathrm{Time}\quad[\mathrm{h}]",
                      ylabel=L"\mathrm{Counts}",dpi=:500,
                      lw=2,c=:black,labelfontsize=14,legendfontsize=10,
		      title="(b)",titlefontsize=10,title_location=:left)
	plot!(t4n,bins=(0:0.5:30),st=:stephist,c=:grey,lw=2,label=L"\tau_{400}")


	lay = @layout [a b]
	p_out = plot(p1,p2,layout=lay,dpi=:300,size=(1000,370),left_margin=5mm,bottom_margin=5mm)

	if save==true
		savefig(p_out,savepath*"Fig4.png")
	else
		return p_out	
	end
	
end

#"Figure 5: P and T at beginning and end of ascent"
function plot_Fig5(save=false)
	t0 = dsnn["t"][:,1] .- 273.15
	p0 = dsnn["p"][:,1] ./ 100

	p_out1 = plot(1:0.5:48,t6_binning_fine(t0,t6n,nm.mean),c=:black,label=false,
		     lw=1.5,ylabel=L"T(0)\quad[°\mathrm{C}]",xlabel=taulabel,labelfontsize=14,
		     fillrange=(t6_prcntl(t0,t6n,10),t6_prcntl(t0,t6n,90)),fillalpha=0.3,
		     title="(a)",titlefontsize=10,title_location=:left)

	plot!(1:0.5:48,t6_binning_fine(t0,t6n,nm.median),c=:black,label=false,lw=1,ls=:dash)
	
	plot!(twinx(),1:0.5:48,t6_binning_fine(p0,t6n,mean),c=:red,label=false,lw=1.5,
	       ls=:solid,ylabel=L"p(0) \quad [\mathrm{hPa}]",
	       guidefontcolor=:red,labelfontsize=15,legendfontsize=10,
	       y_foreground_color_axis=:red,y_foreground_color_text=:red,
	       y_foreground_color_border=:red,ylims=(900,1020),yticks=collect(880:20:1020),
	       size=(500,350),fillrange=(t6_prcntl(p0,t6n,10),t6_prcntl(p0,t6n,90)),fillalpha=0.3)
	plot!(twinx(),1:0.5:48,t6_binning_fine(p0,t6n,median),c=:red,label=false,lw=1,ls=:dash,
	       ylims=(900,1020),ticks=false,axis=false)

	
	t1 = dsnn["t"][:,end] .- 273.15
	p1 = dsnn["p"][:,end] ./ 100

	p_out2 = plot(1:0.5:48,t6_binning_fine(t1,t6n,nm.mean),c=:black,label=false,
		     lw=1.5,ylabel=L"T(1)\quad[°\mathrm{C}]",xlabel=taulabel,labelfontsize=14,
		     fillrange=(t6_prcntl(t1,t6n,10),t6_prcntl(t1,t6n,90)),fillalpha=0.3,
		     ylims=(-60,-28),yticks=collect(-60:5:-30),
		     title="(b)",titlefontsize=10,title_location=:left)
	
	plot!(1:0.5:48,t6_binning_fine(t1,t6n,nm.median),c=:black,label=false,lw=1,ls=:dash)

	plot!(twinx(),1:0.5:48,t6_binning_fine(p1,t6n,mean),c=:red,label=false,lw=2,
	       ls=:solid,ylabel=L"p(1) \quad [\mathrm{hPa}]",
	       guidefontcolor=:red,labelfontsize=15,legendfontsize=10,
	       y_foreground_color_axis=:red,y_foreground_color_text=:red,
	       y_foreground_color_border=:red,ylims=(200,400),yticks=collect(200:25:400),
	       size=(500,350),
	       fillrange=(t6_prcntl(p1,t6n,10),t6_prcntl(p1,t6n,90)),fillalpha=0.3)
	plot!(twinx(),1:0.5:48,t6_binning_fine(p1,t6n,median),c=:red,label=false,lw=1,ls=:dash,
	       ylims=(200,400),ticks=false,axis=false)
        
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1),legend=:top,legend_columns=2)
	
	plot!([],c=:snow4,lw=1.5,label="mean")
	plot!([],c=:snow4,lw=1.5,ls=:dash,label="median")

	lay = @layout [a{0.07h}
			b c]
	p_out = plot(blankp,p_out1,p_out2,layout=lay,dpi=:400, size=(950,350), left_margin=5mm,
		      bottom_margin=5mm,right_margin=5mm)
	
	if save==true
		savefig(p_out,savepath*"Fig5.png")
	else
		return p_out	
	end
end

#"Figure 6: total moisture content at beginning and end of ascent"
function plot_Fig6(save=false)
	qv0 = dsnn["qv"][:,1] .* 1000
	qhy0 = sum([dsnn["q$(hy)"][:,1] for hy in hymet]) .* 1000
	qtot0 = qv0 .+ qhy0
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1))
	
	p1 = plot(1:0.5:48,t6_binning_fine(qtot0,t6n,nm.mean),lw=1,
		     label=false,c=:black,ylims=(6,15),yticks=collect(6:15),
		     ylabel=L"Q_\mathrm{tot/vap}(0) \quad [\mathrm{g/kg}]",xlabel=taulabel,
		     labelfontsize=11,dpi=:400,title="(a)",title_position=:left,
		     fillrange=(t6_prcntl(qtot0,t6n,10),t6_prcntl(qtot0,t6n,90)),fillalpha=0.3,
		     right_margin=5mm,titlefontsize=8)

	plot!(1:0.5:48,t6_binning_fine(qtot0,t6n,nm.median),lw=1,ls=:dash,c=:black,label=false)

	plot!(1:0.5:48,t6_binning_fine(qv0,t6n,mean),lw=1,c=:lawngreen,label=false,
	      fillrange=(t6_prcntl(qv0,t6n,5),t6_prcntl(qv0,t6n,90)),fillalpha=0.3)
	
	plot!(1:0.5:48,t6_binning_fine(qv0,t6n,nm.median),lw=1,ls=:dash,c=:lawngreen,label=false)
	
	plot!(twinx(),1:0.5:48,t6_binning_fine(qhy0,t6n,mean),c=:blue,
	       label=false,lw=1,
               ls=:solid,ylabel=L"Q_\mathrm{hyd}(0) \quad [\mathrm{g/kg}]",
               guidefontcolor=:blue,labelfontsize=11,legendfontsize=11,
               y_foreground_color_axis=:blue,y_foreground_color_text=:blue,
               y_foreground_color_border=:blue,ylims=(-0.005,0.32),yticks=collect(0:0.1:0.5),
	       fillrange=(t6_prcntl(qhy0,t6n,10),t6_prcntl(qhy0,t6n,90)),fillalpha=0.3)

	plot!(twinx(),1:0.5:48,t6_binning_fine(qhy0,t6n,median),lw=1,ls=:dash,
		c=:blue,label=false,ylims=(-0.005,0.32),ticks=false,axis=false)
	
	plot!(blankp,[],[],color=:black,label=L"Q_\mathrm{tot}(0)",lw=1)
	plot!([],[],color=:lawngreen,label=L"Q_\mathrm{vap}(0)",lw=1)
	plot!([],[],color=:blue,label=L"Q_\mathrm{hyd}(0)",lw=1)
	plot!([],[],color=:snow4,label=L"\mathrm{mean}",lw=1)
	plot!([],[],color=:snow4,label=L"\mathrm{median}",lw=1,ls=:dash)
	plot!(legend_columns=5,legendfontsize=9,legend=:top)
	
	qv1 = dsnn["qv"][:,end] .* 1000
	qhy1 = sum([dsnn["q$(hy)"][:,end] for hy in hymet]) .* 1000
	qtot1 = qv1 .+ qhy1

	p2 = plot(1:0.5:48,t6_binning_fine(qtot1,t6n,mean),lw=1,c=:black,label=false,
		     xlabel=taulabel,ylabel=L"Q(1)\quad [\mathrm{g/kg}]",labelfontsize=11,
		     fillrange=(t6_prcntl(qtot1,t6n,10),t6_prcntl(qtot1,t6n,90)),fillalpha=0.3,
		     title="(b)",title_position=:left,titlefontsize=8)

	plot!(1:0.5:48,t6_binning_fine(qtot1,t6n,median),lw=1,c=:black,label=false,ls=:dash)

	plot!(1:0.5:48,t6_binning_fine(qhy1,t6n,mean),lw=1,c=:blue,label=false,
	      fillrange=(t6_prcntl(qhy1,t6n,10),t6_prcntl(qhy1,t6n,90)),fillalpha=0.3)
	plot!(1:0.5:48,t6_binning_fine(qhy1,t6n,median),lw=1,c=:blue,label=false,ls=:dash)

	plot!(1:0.5:48,t6_binning_fine(qv1,t6n,mean),lw=1,c=:lawngreen,label=false,
	      fillrange=(t6_prcntl(qv1,t6n,10),t6_prcntl(qv1,t6n,90)),fillalpha=0.3)
	plot!(1:0.5:48,t6_binning_fine(qv1,t6n,median),lw=1,c=:lawngreen,label=false,ls=:dash)
	
	plot!(dpi=:300)	

	lay = @layout [ a{0.1h}
		       b c ]
	p_out = plot(blankp,p1,p2,layout=lay,dpi=:400,size=(800,350),
		      left_margin=5mm,bottom_margin=5mm)
	
	if save==true
                savefig(p_out,savepath*"Fig6.png")
        else
                return p_out
        end
end

#"Figure 7: RHi histogram at end of ascent and correlation with qv_sati"
#	"Note: uses another function defined directly below this one"
function plot_Fig7(save=false)
	RH_i1= RH_i.(dsnn["t"][:,end],dsnn["qv"][:,end],dsnn["p"][:,end])

	bins1=80:0.5:130

	p1 = plot(RH_i1,st=:stephist,normalize=:true,c=:black,label="all",bins=bins1,
	      xlabel=L"\mathrm{RH_i}(1)\quad [\%]",ylabel="Normalised counts",
	      labelfontsize=14,lw=1.5,right_margin=3mm,left_margin=5mm,bottom_margin=5mm)
	plot!(RH_i1[M_c],st=:stephist,normalize=:true,c=:blue,lw=1.5,
	       label="convective",bins=bins1)
	plot!(RH_i1[M_s],st=:stephist,normalize=:true,c=:coral,lw=1.5,
	       label="slow",bins=bins1)
	plot!(p1,title="(a)",title_location=:left)

	p2 = plot_Qv_vs_Qvsat()
	plot!(p2,title="(b)",title_location=:left,colorbar=false)

	cbar = heatmap((0:0.01:1) .* ones(101,1), legend=:none, xticks=:none,
				  yticks=(1:10:101, string.(0:0.01:0.10)),
				  ylabel="Normalised number of trajectories")

	lay = @layout [ a b c{0.03w}]
	p_out = plot!(p1,p2,cbar,layout=lay,dpi=:300,legendfontsize=10,size=(1000,400))
	if save==true
		savefig(p_out,savepath*"Fig7.png")
	else
		return p_out	
	end
end

#"Used by function above for Figure 7"
function plot_Qv_vs_Qvsat(save=false)
	ar1 = dsnn["qv"][:,end] .* 1000
	ar2 = qv_sat_i.(dsnn["t"][:,end],dsnn["p"][:,end]) .* 1000
	xbins = 0:0.01:1.25
	ybins = 0:0.01:1.25

	xlabel=L"Q_\mathrm{v}(1) \quad [\mathrm{g/kg}]"
	ylabel=L"Q_\mathrm{v,sat,ice}(T(1),p(1)) \quad [\mathrm{g/kg}]"
	loc = :topleft
	
	p_out = plot_hist2d_2ars_general(ar1,ar2,xbins,ybins,xlabel,ylabel,loc,false)
	plot!(xbins,xbins,c=:cyan2,lw=2.2,label="x=y")
	
	return p_out
end


#"Figure 8: plot with CR, PE, DR etc"
#"Note: uses function defined after this one"
function plot_Fig8(save=false)
	p1 = plot_DR()
	plot!(legend=false,cbar=:false,ylims=(0.925,1),yticks=collect(0.92:0.01:1),
	       right_margin=-5mm)
	annotate!(7,0.935,"(a)",annotationcolor=:white,fontsize=:9)
	p2 = plot_DR_mphys()
	plot!(legend=false,cbar=:false,ylims=(0.925,1),yticks=collect(0.92:0.01:1),
	       right_margin=-5mm)
	annotate!(7,0.935,"(b)",annotationcolor=:white,fontsize=:9)
	p3 = plot_DR_mix()
	plot!(legend=false,cbar=:false)
	annotate!(7,-0.02,"(c)",annotationcolor=:white,fontsize=:9)
	p4 = plot_CR()
	plot!(legend=false,cbar=:false,yticks=collect(0.92:0.01:1),ylims=(0.925,1),
	       right_margin=-5mm)
	annotate!(7,0.935,"(d)",annotationcolor=:white,fontsize=:9)
	p5 = plot_PE()
	plot!(legend=false,cbar=:false,right_margin=-5mm)
	annotate!(7,0.987,"(e)",annotationcolor=:white,fontsize=:9)
	p6 = plot_P1_over_qtot0()
	plot!(legend=false,cbar=:false,yticks=collect(0.4:0.2:1.2))
	annotate!(7,0.49,"(f)",annotationcolor=:white,fontsize=:9)
	cbar = heatmap((0:0.01:1)' .* ones(101,1), legend=:none, yticks=:none,
				   xticks=(1:10:101, string.(0:0.01:0.1)),
				   xlabel="Normalised number of trajectories")
	lay = @layout [a b c 
			 d e f
			 g{0.05h} ]
	p_out = plot(p1,p2,p3,p4,p5,p6,cbar,layout=lay,dpi=:400,size=(1000,600),
		      left_margin=2.5mm,bottom_margin=3mm,right_margin=-3mm)
	if save==false
		return p_out
	elseif save==true	
		savefig(p_out,savepath*"Fig8.png")
	end
end

#"Figure 8 uses the following functions:"
#-------------------------------------------
function plot_PE()
	ar = dsPE["PE"][:,end]
	ybins = 0.985:0.0003:1
	ylabel = L"\mathrm{PE}"
	loc = :bottom
	h_ar = histogram2d((t6n,ar),clims=(0,0.1),xlims=(0,48),ylims=(ybins[1],ybins[end]),
	                        xlabel=taulabel,ylabel=ylabel,
                                bins=(0:1:48,ybins),labelfontsize=14,
                                colorbartitle="\n Normalised number of trajectories",
                                right_margin=5mm,left_margin=2mm,
                                legend=loc,c=cgrad(:inferno))
	plot!(1:0.5:48,t6_binning_fine(ar,t6n,nm.mean),label="mean",lw=1.5,c=:cyan)
	plot!(1:0.5:48,t6_binning_fine(ar,t6n,nm.median),label="median",ls=:dash,lw=1.5,c=:cyan)
	plot!(h_ar,yticks=collect(0.985:0.005:1))

	normalize_hist2d(h_ar)

	return h_ar
end

function plot_DR()
	ar = dsPE["DR"][:,end]
	ybins = 0.925:0.001:1
	ylabel = L"\mathrm{DR}"
	loc=:bottomleft
	p_out = plot_hist2d_t6_general(ar,ybins,ylabel,loc)
	plot!(yticks=[0.94,0.96,0.98,1],dpi=:300)
	return p_out
end

function plot_DR_mphys()
	ar = dsPE["DR_mphys"][:,end]
	ybins = 0.925:0.0021:1
	ylabel=L"\mathrm{DR_{mphys}}"
	loc = :bottom
	p_out = plot_hist2d_t6_general(ar,ybins,ylabel,loc)
	plot!(dpi=:300)
	return p_out
end

function plot_DR_mix()
	ar = dsPE["DR_mix"][:,end]
	ybins = -0.1:0.015:0.5
	ylabel=L"\mathrm{DR_{mix}}"
	loc = :topleft
	h_ar = plot_hist2d_t6_general(ar,ybins,ylabel,loc)
        return h_ar
end


function plot_CR()
	ar = dsPE["CR"][:,end]
	ybins=0.925:0.0008:1
	ylabel=L"\mathrm{CR}"
	loc = :bottomleft
	p_out = plot_hist2d_t6_general(ar,ybins,ylabel,loc)
	return p_out
end

function plot_P1_over_qtot0()
	qtot0 = dsPE["qhy"][:,1] .+ dsnn["qv"][:,1]
	P1 = dsPE["allfl"][:,end]
	ar = P1 ./ qtot0

	ybins=0.4:0.02:1.2
	ylabel=L"P(1)/Q_\mathrm{tot}(0)"
	loc=:topright
	p_out = plot_hist2d_t6_general(ar,ybins,ylabel,loc)
	return p_out
end
#-------------------------------------------

#"Figure 9: CR, Temperature, over normed ascent time, and RH etc in temperature bins"
function plot_Fig9(save=false)
	p1 = plot_twotau_general(Array(dsPE["CR"]))
	plot!(p1,ylabel=L"\mathrm{CR}",legend=false,title="(a)",xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")
	
	p2 = plot_twotau_general(Array(dsnn["t"]) .- 273.15)
	plot!(p2,ylabel=L"\mathrm{Temperature} \quad [°\mathrm{C}]",legend=false,title="(b)",ylims=(-60,22),
		   xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")

	RHw = RH_w.(Array(dsnn["t"]), Array(dsnn["qv"]), Array(dsnn["p"]))
	p3 = plot_Tbinned_var_twotau(RHw)
	plot!(p3,ylabel=L"\mathrm{RH}_\mathrm{w} \quad [\%]",xlims=(-40,20),legend=false,title="(c)")
	
	GC.gc()
	
	RHi = RH_i.(Array(dsnn["t"]), Array(dsnn["qv"]), Array(dsnn["p"]))
	p4 = plot_Tbinned_var_twotau(RHi)
	plot!(p4,ylabel=L"\mathrm{RH}_\mathrm{i} \quad [\%]",xlims=(-50,0),legend=false,title="(d)")
	
	GC.gc()

	liq_frac = (Array(dsnn["qc"]) .+ Array(dsnn["qr"])) ./ Array(dsPE["qhy"])

	p5 = plot_Tbinned_var_twotau(liq_frac)
	plot!(p5,ylabel=L"\mathrm{Liquid\,\,\,fraction}",legend=false,ylims=(0,1),xlims=(-40,20),title="(e)")
	
	p6 = plot_Tbinned_var_twotau(Array(dsnn["qni"]), true)
	plot!(p6,legend=false,ylabel=L"\mathrm{qni} \quad [\mathrm{1/kg}]",yscale=:log10,ylims=(10,1e7),title="(f)")
	
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
			   yticks=false,xlim=(0,1),ylim=(0,1))
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=2)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=2,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=6)
	
	lay = @layout [a{0.09h}
		       b c
		       d e
		       f g]
	p_out = plot(blankp,p1,p2,p3,p4,p5,p6,layout=lay,
		      dpi=:400,labelfontsize=8,tickfontsize=5,titlefontsize=6,
		      title_position=:left,size=(650,500),lw=1)

	if save==true
                savefig(p_out,savepath*"Fig9.png")
        else
                return p_out
        end
end	

#"Figure 10: microphysics"
function plot_Fig10(save=false)
	p1 = plot_twotau_general(Array(dsPE["PE"]))
	plot!(p1,ylabel=L"\mathrm{PE}",legend=false,title="(a)",xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")
	
	p2 = plot_twotau_general(Array(dsPE["allfl"]) .* 1000)
	plot!(p2,ylabel=L"\mathrm{Rain \,\,\, precip.}\quad [\mathrm{g/kg}]",legend=false,title="(b)",
		   xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")

	GC.gc()

	p3 = plot_twotau_general(Array(dsPE["coldf"]) .* 1000)
	plot!(p3,ylabel=L"\mathrm{Frozen \,\,\, precip.}\quad [\mathrm{g/kg}]",legend=false,title="(c)",
		   xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")

	p4 = plot_twotau_general(Array(dsPE["qhy"]) .* 1000)
	plot!(p4,ylabel=L"\mathrm{Total \,\,\, hydr.}\quad [\mathrm{g/kg}]",legend=false,title="(d)",
		   xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")

	GC.gc()
	
	p5 = plot_twotau_general(Array(dsnn["qs"]) .* 1000)
	plot!(p5,ylabel=L"\mathrm{qs}\quad [\mathrm{g/kg}]",legend=false,title="(e)",
		   xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}",ylims=(0,0.7),yticks=collect(0:0.1:0.7))

	p6 = plot_twotau_general(Array(dsnn["qg"]) .* 1000)
	plot!(p6,ylabel=L"\mathrm{qg}\quad [\mathrm{g/kg}]",legend=false,title="(f)",
		   xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}",ylims=(0,0.7),yticks=collect(0:0.1:0.7))

	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
			   yticks=false,xlim=(0,1),ylim=(0,1))
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=2)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=2,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=6)
	
	lay = @layout [a{0.09h}
			   b c d
			   e f g]
	p_out = plot(blankp,p1,p2,p3,p4,p5,p6,layout=lay,
			  dpi=:400,labelfontsize=8,tickfontsize=5,titlefontsize=6,
			  title_position=:left,size=(700,450),lw=1)
	
	if save==true
				savefig(p_out,savepath*"Fig10.png")
		else
				return p_out
		end
end

end

#"Figure 11: qni qi ri RH after asc" 
function plot_Fig11(save=false)
	qni_all = Dict()
	qi_all = Dict()
	ri_all = Dict()
	alt_all = Dict()
	nihh_all = Dict()
	nihom_all = Dict()
	altr = Array(dsnr["alt"])
	tr = Array(dsnr["t"])
	qvr = Array(dsnr["qv"])
	pr = Array(dsnr["p"])
	qnir = Array(dsnr["qni"])
	qir = Array(dsnr["qi"])
	RHr = RH_i.(tr,qvr,pr)
	rir = calc_radius_ice(qir,qnir)
	GC.gc()

	RH_all = Dict()
        ixd = Dict()
        mxd = Dict()
        xar = 0:1:45

	for i in xar
                print("$(i), ")
                ix,mx = get_tWCB_plus_ts_indices(i)
                ixd[i] = ix
                mxd[i] = mx
		alt_all[i] = altr[ix[mx]]
		
		qni_all[i] = qnir[ix[mx]]	
		qi_all[i] = qir[ix[mx]] .* 1000
		ri_all[i] = rir[ix[mx]] .* 1000
		msk = (qi_all[i] .< 1e-10) .| (qni_all[i] .< 5)
		ri_all[i][msk] .= NaN
		RH_all[i] = RHr[ix[mx]]
		GC.gc()

	end
	
	p1 = plot_after_asc(xar,qni_all,mxd)
	plot!(ylabel="qni [1/kg]",legend=false,ylims=(10,1e7),yticks=[1e1,1e2,1e3,1e4,1e5,1e6,1e7],
	       yscale=:log10,title="(a)",title_location=:left,titlefontsize=9,labelfontsize=9)
	
	p2 = plot_after_asc(xar,qi_all,mxd)
	plot!(ylabel="qi [g/kg]",legend=false,ylims=(1e-6,1),yticks=[1e-6,1e-4,1e-2,1],
	       yscale=:log10,title="(b)",title_location=:left,titlefontsize=9,labelfontsize=9)
	
	p3 = plot_after_asc(xar,ri_all,mxd)
	plot!(ylabel="Mean geometric radius [mm]",legend=false,ylims=(0.0,0.2),yticks=collect(0.0:0.05:0.2),
	       title="(c)",title_location=:left,titlefontsize=9,labelfontsize=9)
	
	p4 = plot_after_asc(xar,RH_all,mxd)
	plot!(ylabel="RH over ice [%]",legend=false,ylims=(40,120),yticks=collect(40:20:120),
	       title="(d)",title_location=:left,titlefontsize=9,labelfontsize=9)
	
	
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1))
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=1.5)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=1.5,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=9)
	
	lay = @layout [a{0.1h}
			b c
			 d e]
	p_out = plot(blankp,p1,p2,p3,p4,layout=lay,dpi=:300,size=(800,600),left_margin=2mm)
        
	if save==true
                savefig(p_out,savepath*"Fig11.png")
        else
                return p_out
        end
end


#----------------------------------------------------------------------------
#	"Appendix"
#----------------------------------------------------------------------------

#"Figure C1: dp2 plus example trajectories"
function plot_FigC1(save=false)
	ar1 = t6n
	dp2 = find_max_dp2()
	ar2 = dp2
	xbins = 0:0.5:48
	ybins = 50:5:800
        xlabel=taulabel
	ylabel=L"\mathrm{max}(\Delta p_\mathrm{2\,h})\quad [\mathrm{hPa}]"
	loc = :topright
	p_out1 = plot_hist2d_2ars_general(ar1,ar2,xbins,ybins,xlabel,ylabel,loc,true)
	#plot!([8,48,48,25,15,10,8,8],[600,600,300,300,350,450,550,600],label=false,lw=2,c=:coral)
	plot!([20,48,48,20,20],[350,350,600,600,350],label=false,lw=2,c=:coral)
	annotate!(35,570,text("Excluded",15,color=:coral))
	title!("(a)",titlefontsize=:10,title_location=:left)
	plot!(colorbar=false,right_margin=-3mm)
	
	cbar = heatmap((0:0.01:1)' .* ones(101,1), legend=:none, yticks=:none,
				   xticks=(1:10:101, string.(0:0.01:0.1)),
				   xlabel="Normalised number of trajectories")
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1),legend=:false,legend_columns=2)
	
	k1 = 924 #idx = 4218
	k2 = 341 #318136 #idx = 1861466

	i1 = find_start(k1):find_stop(k1)
	i2 = find_start(k2):find_stop(k2)
	j1 = (an[k1] + ns[k1]):(an[k1] + ns[k1] + tn[k1])
	j2 = (an[k2] + ns[k2]):(an[k2] + ns[k2] + tn[k2])
	rt = Array(dsnr["rtime"])
	pnr1 = dsnr["p"][k1,:] ./ 100
	pnr2 = dsnr["p"][k2,:] ./ 100 
	p_out2 = plot(rt,pnr1,yflip=:true,lw=2,c=:black,
		      label=L"\mathrm{total}",
		      xlabel=L"\mathrm{Sim.Time}\quad[\mathrm{h}]",
		      legend=:topleft,legendfontsize=10,legendcolumn=1,
		      ylabel=L"\mathrm{Pressure}\quad[\mathrm{hPa}]",
		      labelfontsize=14)
	plot!(rt[i1],pnr1[i1],c=:orange,lw=2,label=L"\tau_\mathrm{WCB}")
	plot!(rt[j1],pnr1[j1],c=:red,lw=2,label=L"\tau_{600}")
	plot!(rt,pnr2,yflip=:true,lw=2,c=:black,label=false)
	plot!(rt[i2],pnr2[i2],c=:orange,lw=2,label=false)
	plot!(rt[j2],pnr2[j2],c=:red,lw=2,label=false)
	title!("(b)",titlefontsize=:10,title_location=:left)

	lay = @layout [grid(2,1,heights=[0.96,0.04]) b]

	p_out = plot(p_out1,cbar,p_out2,layout=lay,dpi=:300,left_margin=5mm,bottom_margin=5mm,
		      size=(1000,500))
	
	if save==true
		savefig(p_out,savepath*"FigC1.png")
	else
		return p_out	
	end
end

#"Figure C2: wind shear"
#"The data in the arrays seen here is obtained by running calculate_mean_wind_shear.jl"
function plot_FigC2(save=false)
	EW = [0.382, 0.383, 0.376, 0.388, 0.408, 0.423, 0.443, 0.463, 0.474, 0.486, 
	     0.5, 0.51, 0.525, 0.535, 0.555, 0.559, 0.573, 0.589, 0.596, 0.601, 0.609, 
	     0.619, 0.627, 0.621, 0.636, 0.636, 0.648, 0.657, 0.657, 0.667, 0.67, 0.679, 
	     0.664, 0.67, 0.675, 0.655, 0.653, 0.652, 0.659, 0.655, 0.64, 0.651, 0.64, 
	     0.647, 0.666, 0.673, 0.654, 0.662, 0.681, 0.65, 0.659, 0.641, 0.647, 0.656, 
	     0.644, 0.633, 0.626, 0.622, 0.61, 0.616, 0.592, 0.608, 0.574, 0.58, 0.568, 
	     0.591, 0.558, 0.572, 0.564, 0.531, 0.538, 0.517, 0.519, 0.529, 0.51, 0.525, 
	     0.519, 0.518, 0.501, 0.5, 0.51, 0.499, 0.481, 0.488, 0.477, 0.47, 0.471, 
	     0.475, 0.468, 0.463, 0.458, 0.453, 0.461, 0.442, 0.433]
	WZ = [0.006578, 0.007236, 0.008566, 0.009135, 0.010108, 0.010532, 0.010805, 
	      0.01097, 0.010794, 0.011068, 0.010843, 0.010919, 0.011009, 0.010893, 0.010763, 
	      0.010941, 0.010871, 0.011271, 0.011055, 0.010964, 0.01085, 0.010759, 0.010757, 
	      0.01044, 0.01039, 0.010382, 0.01052, 0.010407, 0.01047, 0.010443, 0.010588, 
	      0.010464, 0.010168, 0.010212, 0.009965, 0.009814, 0.009824, 0.009578, 0.00982, 
	      0.009458, 0.009607, 0.009623, 0.009534, 0.009562, 0.009555, 0.009654, 0.009509, 
	      0.009457, 0.009769, 0.009298, 0.009141, 0.009373, 0.009371, 0.00928, 0.009511, 
	      0.009315, 0.009427, 0.00888, 0.009187, 0.008752, 0.00884, 0.008739, 0.008659, 
	      0.008569, 0.008895, 0.008841, 0.008648, 0.008458, 0.009079, 0.008351, 0.008543, 
	      0.008343, 0.008301, 0.008509, 0.008078, 0.008549, 0.008226, 0.00829, 0.008049, 
	      0.007969, 0.007631, 0.008146, 0.007412, 0.007988, 0.007627, 0.007601, 0.007485, 
	      0.007345, 0.007477, 0.00721, 0.006984, 0.007438, 0.007425, 0.006492, 0.006975] 
	
	p_out = plot(1:0.5:48,EW,xlabel=taulabel,ylabel="Mean E-W wind shear [1/s]",lw=1.5,
		  c=:black,label=false)
	plot!(twinx(),1:0.5:48,WZ,ylabel="Mean vertical wind shear [1/s]",lw=1.5,c=:red,
	       label=false,guidefontcolor=:red,y_foreground_color_axis=:red,
	       y_foreground_color_text=:red,y_foreground_color_border=:red)
	plot!(labelfontsize=:10,dpi=:300)

	if save==true
                savefig(p_out,savepath*"FigC2.png")
        else
                return p_out
        end
end

#"Figure C3: precipitation types"
function plot_FigC3(save=false)
	P1 = sum([flux_fast_last(hy,"normed") for hy in hymet[2:end]]) .* 1000
	Pr = flux_fast_last("r","normed") .* 1000
	Pg = flux_fast_last("g","normed") .* 1000
	Ps = flux_fast_last("s","normed") .* 1000
	Pice = flux_fast_last("i","normed") .* 1000
	Ph = flux_fast_last("h","normed") .* 1000

	Pcold = Pg .+ Ps .+ Pice .+ Ph

	p1 = plot(1:0.5:48,t6_binning_fine(P1,t6n,mean),lw=1.5,c=:black,
		     label=L"P_\mathrm{total}(1)",xlabel=taulabel,ylabel="Flux [g/kg]",
		     legend=:outertop,legendfontsize=:9,legend_columns=3,labelfontsize=12,
		     fillrange=(t6_prcntl(P1,t6n,10),t6_prcntl(P1,t6n,90)),fillalpha=0.3)
	plot!(1:0.5:48,t6_binning_fine(Pr,t6n,mean),lw=1.5,ls=:solid,c=:blue,
	      label=L"P_\mathrm{rain}(1)",
	      fillrange=(t6_prcntl(Pr,t6n,10),t6_prcntl(Pr,t6n,90)),fillalpha=0.3)
	plot!(1:0.5:48,t6_binning_fine(Pcold,t6n,mean),lw=1.5,ls=:solid,c=:turquoise2,
	      label=L"P_\mathrm{frozen}(1)",
	      fillrange=(t6_prcntl(Pcold,t6n,10),t6_prcntl(Pcold,t6n,90)),fillalpha=0.2)
	hline!([0],label=false,c=:black)
	plot!(title="(a)",title_location=:left,titlefontsize=8)
	
	p2 = plot(1:0.5:48,t6_binning_fine(Pice,t6n,mean),lw=1.5,ls=:solid,c=:black,
	      label=L"P_\mathrm{ice}(1)",legend=:outertop,legendfontsize=:9,
	      legend_columns=4,labelfontsize=12,xlabel=taulabel,ylabel="Frozen flux [g/kg]",
	      fillrange=(t6_prcntl(Pice,t6n,10),t6_prcntl(Pice,t6n,90)),fillalpha=0.5)
	hline!([0],label=false,c=:black)
	plot!(title="(b)",title_location=:left,titlefontsize=8)
	
	plot!(1:0.5:48,t6_binning_fine(Pg,t6n,mean),lw=1.5,ls=:solid,c=:orangered,
	      label=L"P_\mathrm{graupel}(1)",
	      fillrange=(t6_prcntl(Pg,t6n,10),t6_prcntl(Pg,t6n,90)),fillalpha=0.4)
	
	plot!(1:0.5:48,t6_binning_fine(Ps,t6n,mean),lw=1.5,ls=:solid,c=:cyan,
	      label=L"P_\mathrm{snow}(1)",
	      fillrange=(t6_prcntl(Ps,t6n,10),t6_prcntl(Ps,t6n,90)),fillalpha=0.3)

	plot!(1:0.5:48,t6_binning_fine(Ph,t6n,mean),lw=1.5,ls=:dash,c=:yellow,
	      label=L"P_\mathrm{hail}(1)",
	      fillrange=(t6_prcntl(Ph,t6n,10),t6_prcntl(Ph,t6n,90)),fillalpha=0.5)

	
	lay = @layout [ a b ]
	p_out = plot(p1,p2,layout=lay,dpi=:300,size=(800,400),left_margin=4mm,bottom_margin=4mm,
		      right_margin=5mm)
	if save==true
                savefig(p_out,savepath*"FigC3.png")
        else
                return p_out
        end
end

#"Figure C4: DR mix contributions"
function plot_FigC4(save=false)
	DR_mix = dsPE["DR_mix"][:,end]
	qtot0 = dsPE["qhy"][:,1] .+ dsnn["qv"][:,1]
	
	turbs = sum([dsnn["q$(hy)turb"][:,end] for hy in ["c","v"]])
	Tu = -turbs ./ qtot0
	
	convs = sum([dsnn["q$(x)conv"][:,end] for x in ["v","c","r","i","s"]])
	Co = -convs ./ qtot0

	res = calc_res(101)
	R = res ./ qtot0

	ylims=(0,0.5)
	p_out = plot(1:0.5:48,t6_binning_fine(DR_mix,t6n,mean),lw=1.5,c=:black,
		     label=L"\mathrm{DR_{mix}}",xlabel=taulabel,dpi=:300,
		     ylabel="Fractional removal of moisture",
		     legend_position=(0.2,0.96),legend_columns=1,ylims=ylims,
		     legendfontsize=:10,legend_background_color=nothing,
		     foreground_color_legend=nothing,
		     fillrange=(t6_prcntl(DR_mix,t6n,10),t6_prcntl(DR_mix,t6n,90)),
		     fillalpha=0.2)
	
	plot!(1:0.5:48,t6_binning_fine(Tu,t6n,mean),lw=2,c=:royalblue1,ls=:solid,
	      label=L"Q_\mathrm{turb}(1)/Q_\mathrm{tot}(0)",
	      fillrange=(t6_prcntl(Tu,t6n,10),t6_prcntl(Tu,t6n,90)),fillalpha=0.5)

	plot!(1:0.5:48,t6_binning_fine(Co,t6n,mean),lw=2,c=:magenta,ls=:solid,
	      label=L"Q_\mathrm{conv}(1)/Q_\mathrm{tot}(0)",
	      fillrange=(t6_prcntl(Co,t6n,10),t6_prcntl(Co,t6n,90)),fillalpha=0.5)

	plot!(1:0.5:48,t6_binning_fine(R,t6n,mean),lw=2,c=:orange,ls=:solid,
	      label=L"\mathfrak{R}(1)/Q_\mathrm{tot}(0)",
	      fillrange=(t6_prcntl(R,t6n,10),t6_prcntl(R,t6n,90)),fillalpha=0.4)
 
	if save==true
		savefig(p_out,savepath*"FigC4.png")
	else
		return p_out	
	end
end


#"Figure C5: DR_mix, pressure and vertical velocity over time"
function plot_FigC5(save=false)
	p1 = plot_twotau_general(Array(dsPE["DR_mix"]))
	plot!(p1,ylabel=L"\mathrm{DR}_\mathrm{mix}",legend=false,title="(a)",xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")
		
	p2 = plot_twotau_general(Array(dsnn["p"]) ./ 100)
		plot!(p2,ylabel=L"\mathrm{Pressure}\quad [\mathrm{hPa}]",legend=false,ylims=(200,1000),
			   title="(b)",yflip=true,xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")
	
	wv = zeros(size(Array(dsPE["qhy"])))
    wv[:,2:end] = calc_vertical_velocity()

    p3 = plot_twotau_general(wv)
    plot!(p3,ylabel=L"W_v \quad [\mathrm{m/s}]",legend=false,ylims=(0,1.4),
               title="(c)",yticks=collect(0:0.2:1.4),xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")

	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1))
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=1.5)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=1.5,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=6)
	
	lay = @layout [a{0.1h}
			b c d]
	p_out = plot(blankp,p1,p2,p3,layout=lay,dpi=:300,size=(700,300),left_margin=2mm,title_position=:left,
		      labelfontsize=8,tickfontsize=5,titlefontsize=6,legendfontsize=6)

	if save==true
		savefig(p_out,savepath*"FigC5.png")
	else
		return p_out	
	end

end

#"Figure C6: ice processes"
function plot_FigC6(save=false)
	p1 = plot_twotau_general(-(Array(dsnn["evap"]) .- Array(dsnn["r_evap"])) .* 1000)
	plot!(p1,ylabel=L"\mathrm{evap}+\mathrm{rain\,\,evap}\quad[\mathrm{g/kg}]",ylims=(0,3),
		   legend=false,title="(a)",xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")
		
	p2 = plot_twotau_general(Array(dsnn["wbf"]) .* 1000)
		plot!(p2,ylabel=L"\mathrm{WBF}\quad [\mathrm{g/kg}]",legend=false,ylims=(0,3),
			   title="(b)",xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")
	
	p3 = plot_twotau_general(Array(dsnn["depo"]) .* 1000)
		plot!(p3,ylabel=L"\mathrm{depo}\quad [\mathrm{g/kg}]",legend=false,ylims=(0,6),
			   title="(c)",xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")
	
	p4 = plot_twotau_general(Array(dsnn["qx_rim"]) .* 1000)
		plot!(p4,ylabel=L"\mathrm{qx\_rim}\quad [\mathrm{g/kg}]",legend=false,ylims=(0,6),
			   title="(d)",xlabel=L"\mathrm{Normalised \,\,\, ascent \,\,\, time}")

	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
			   yticks=false,xlim=(0,1),ylim=(0,1))
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=1.5)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=1.5,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=6)
	
	lay = @layout [a{0.1h}
			b c
			d e]
	p_out = plot(blankp,p1,p2,p3,p4,layout=lay,dpi=:300,size=(500,400),left_margin=2mm,title_position=:left,
		      labelfontsize=8,tickfontsize=5,titlefontsize=6,legendfontsize=6)

	if save==true
		savefig(p_out,savepath*"FigC6.png")
	else
		return p_out	
	end

end

#-------------------------------------------------------------------------
#	"Plots that identify findings not shown in the Paper"
#-------------------------------------------------------------------------

#"Figure Not shown: Pressure, deposition, potential temperature and PV after the ascent"
function plot_after_asc_appendix(save=false)
	pres_all = Dict()
	depo_all = Dict()
	theta_all = Dict()
	pv_all = Dict()
        ixd = Dict()
        mxd = Dict()
	i0,m0 = get_tWCB_plus_ts_indices(0)
	depo1 = Array(dsnr["depo"])[i0]
		xar = 0:1:45
	for i in xar
				print("$(i), ")
				ix,mx = get_tWCB_plus_ts_indices(i)
				ixd[i] = ix
				mxd[i] = mx
		temp = Array(dsnr["t"])[ix[mx]]
		pres_all[i] = Array(dsnr["p"])[ix[mx]] ./ 100	
		depo_all[i] = (Array(dsnr["depo"])[ix[mx]] .- depo1[mx]) .* 1000
		theta_all[i] = Theta.(temp,pres_all[i] .* 100)
		pv_all[i] = Array(dsnr["pv"])[ix[mx]] .* 1e6
		GC.gc()

	end
	
	p1 = plot_after_asc(xar,pres_all,mxd)
	plot!(ylabel="Pressure [hPa]",legend=false,ylims=(200,450),yticks=collect(200:50:450),
	       yflip=:true,title="(a)",title_location=:left,titlefontsize=9,labelfontsize=9)
	p2 = plot_after_asc(xar,depo_all,mxd)
	plot!(ylabel="depo [g/kg]",legend=false,
	       title="(b)",title_location=:left,titlefontsize=9,labelfontsize=9)
	p3 = plot_after_asc(xar,theta_all,mxd)
	plot!(ylabel="Pot. temp [K]",legend=false,
	       title="(c)",title_location=:left,titlefontsize=9,labelfontsize=9)
	p4 = plot_after_asc(xar,pv_all,mxd)
	plot!(ylabel="PV [PVU]",legend=false,
	       title="(d)",title_location=:left,titlefontsize=9,labelfontsize=9)
	
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1))
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=1.5)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=1.5,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=8)
	
	lay = @layout [a{0.1h}
			b c
			 d e]
	p_out = plot(blankp,p1,p2,p3,p4,layout=lay,dpi=:300,size=(800,600),left_margin=2mm)
        
	if save==true
                savefig(p_out,savepath*"appendix_after_asc.png")
        else
                return p_out
        end

end

#"frozen fraction and SLW during normalised ascent"
function plot_frozen_fraction(save=false)
	q_liq = Array(dsnn["qc"]) .+ Array(dsnn["qr"])
	qhy = Array(dsPE["qhy"])
	temp = Array(dsnn["t"]) .- 273.15
	t0msk = temp .> 0
	fr_frac = 1 .- (q_liq ./ qhy)
	SLW = q_liq .* 1000
	SLW[t0msk] .= NaN
	SLW[:,1:5] .= NaN
	GC.gc()
	p1 = plot_twotau_general(fr_frac)
	p2 = plot_twotau_general(SLW)
	
	plot!(p1,ylabel="Frozen fraction",legend=false,title="(a)",
	       labelfontsize=9,titlefontsize=9,title_location=:left)
	plot!(p2,ylabel="SLW [g/kg]",legend=false,title="(b)",
	       labelfontsize=9,titlefontsize=9,title_location=:left)
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1))
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=1.5)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=1.5,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=8)
	
	lay = @layout [a{0.12h}
			b c]
	p_out = plot!(blankp,p1,p2,layout=lay,dpi=:300,size=(800,350),
		       left_margin=3mm,bottom_margin=3mm)

	if save==true
                savefig(p_out,savepath*"SLW_frozen_fraction.png")
        else
                return p_out
        end

end

#"moisture in pressure bins"
function plot_pbin_moisture(save=false)
	qv = Array(dsnn["qv"]) .* 1000
	qhy = Array(dsPE["qhy"]) .* 1000
	qtot = qv .+ qhy

	GC.gc()
	w_v = Array(dsnn["w_v"])

	p1 = plot_Pbinned_var_twotau(qtot)
	GC.gc()
	p2 = plot_Pbinned_var_twotau(qv)
	GC.gc()
	p3 = plot_Pbinned_var_twotau(qhy)
	GC.gc()
	p4 = plot_Pbinned_var_twotau(w_v)
	GC.gc()

	yt1 = vcat(0.1:0.1:0.5,1:1:4)
	yt2 = vcat(0.1:0.1:0.5,1:1:4)
	yt3 = [1e-3,5e-3,1e-2,5e-2,1e-1,5e-1,1]

	plot!(p1,legend=false,xlims=(250,500),ylims=(0.05,4.5),yscale=:log10,yticks=(yt1,string.(yt1)),
	       ylabel=L"Q_\mathrm{tot}\quad[\mathrm{g/kg}]",title="(a)",left_margin=3mm)
	plot!(p2,legend=false,xlims=(250,500),ylims=(0.05,4.5),yscale=:log10,yticks=(yt2,string.(yt2)),
	       ylabel=L"Q_\mathrm{vap}\quad[\mathrm{g/kg}]",title="(b)",left_margin=3mm)
	plot!(p3,legend=false,xlims=(250,500),ylims=(0.0005,2),yscale=:log10,yticks=(yt3,string.(yt3)),
	       ylabel=L"Q_\mathrm{hyd}\quad[\mathrm{g/kg}]",title="(c)")
	plot!(p4,legend=false,xlims=(250,500),ylims=(0,1),
	       ylabel=L"\mathrm{w_v}\quad[\mathrm{m/s}]",title="(d)")

	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1))
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=1.5)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=1.5,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=8)

	lay = @layout [a{0.1h}
		       b c
		       d e]

	p_out = plot(blankp,p1,p2,p3,p4,layout=lay,labelfontsize=12,title_position=:left,
		      titlefontsize=7,size=(900,600),dpi=:350,left_margin=3mm)
	
	if save==true
                savefig(p_out,savepath*"pbin_moisture.png")
        else
                return p_out
        end
end

function plot_micro_time_new(save=false)
	te = Array(dsnn["t"])
	print("1")
	qv = Array(dsnn["qv"])
	print("2")
	p = Array(dsnn["p"])
	print("3")

	RH_iall = RH_i.(te,qv,p);
	RH_wall = RH_w.(te,qv,p);
	print("4")
	temp = te .- 273.15;
	ts0m = temp .<= 0;
	tb0m = temp .> 0;

	te = nothing
	qv = nothing
	p = nothing

	print("5")
	RH_iall[tb0m] .= 0;
	RH_wall[ts0m] .= 0;
	RH_tot = RH_iall .+ RH_wall;
	print("6")
	allfl = Array(dsPE["allfl"]);
	coldf = Array(dsPE["coldf"]);
	rainf = Array(dsPE["rainf"]);

	qc = Array(dsnn["qc"])
	qr = Array(dsnn["qr"])
	qhy = Array(dsPE["qhy"])
	liq_frac = (qc .+ qr) ./ qhy

	qc = nothing
	qr = nothing
	qhy = nothing

	print("7")
	GC.gc()


	p1 = plot_twotau_general(allfl .* 1000)
	plot!(p1,ylabel="Total prec. [g/kg]",legend=false,ylims=(-1,14),
	       yticks=collect(0:2:14),title="(a)")
		
	p2 = plot_twotau_general(rainf .* 1000)
	plot!(p2,ylabel="Rain prec. [g/kg]",legend=false,ylims=(-1,14),
	       yticks=collect(0:2:14),title="(b)")
		
	p3 = plot_twotau_general(coldf .* 1000)
	plot!(p3,ylabel="Frozen prec. [g/kg]",legend=false,title="(c)")
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1))

	p4 = plot_twotau_general(liq_frac)
	plot!(p4,ylabel="liquid fraction [g/kg]",ylims=(0,1),legend=false,yticks=collect(0:0.2:1),
	       left_margin=2mm,title="(d)")

	p5 = plot_twotau_general(RH_iall)
	plot!(p5,ylabel="RHi [%]",ylims=(75,120),legend=false,yticks=collect(80:10:120),
	       title="(e)")
	
	p6 = plot_twotau_general(RH_wall)
	plot!(p6,ylabel="RHw [%]",ylims=(75,120),legend=false,yticks=collect(80:10:120),
	       title="(f)")
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=2)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=2,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=7)
	GC.gc()
	
	lay = @layout [a{0.09h}
		       b c d
		       e f g]
	p_out = plot(blankp,p1,p2,p3,p4,p5,p6,layout=lay,
		      dpi=:400,labelfontsize=6,tickfontsize=5,titlefontsize=6,
		      title_position=:left,size=(500,400),left_margin=5mm,bottom_margin=5mm)

	if save==true
                savefig(p_out,savepath*"micro_time.png")
        else
                return p_out
        end
end

function plot_micro_time(save=false)
	evap = Array(dsnn["evap"]) .* 1000;
	revap = Array(dsnn["r_evap"]) .* 1000;
	evaps = evap .- revap;
	cnd = Array(dsnn["cond"]) .* 1000;
	depo = Array(dsnn["depo"]) .* 1000;
	qx_rim = Array(dsnn["qx_rim"]) .* 1000;
	RH_iall = RH_i.(Array(dsnn["t"]), Array(dsnn["qv"]), Array(dsnn["p"]));
	RH_wall = RH_w.(Array(dsnn["t"]), Array(dsnn["qv"]), Array(dsnn["p"]));
	temp = Array(dsnn["t"]) .- 273.15;
	ts0m = temp .<= 0;
	tb0m = temp .> 0;
	RH_iall[tb0m] .= 0;
	RH_wall[ts0m] .= 0;
	RH_tot = RH_iall .+ RH_wall;
	allfl = Array(dsPE["allfl"]);
	coldf = Array(dsPE["coldf"]);
	rainf = Array(dsPE["rainf"]);

	wbf = Array(dsnn["wbf"]) .* 1000;

	p1 = plot_twotau_general(allfl .* 1000)
	plot!(p1,ylabel="Total prec. [g/kg]",legend=false,ylims=(-1,14),
		   yticks=collect(0:2:14),title="(a)")
		
	p2 = plot_twotau_general(rainf .* 1000)
	plot!(p2,ylabel="Rain prec. [g/kg]",legend=false,ylims=(-1,14),
		   yticks=collect(0:2:14),title="(b)")
		
	p3 = plot_twotau_general(coldf .* 1000)
	plot!(p3,ylabel="Frozen prec. [g/kg]",legend=false,title="(c)")
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
			   yticks=false,xlim=(0,1),ylim=(0,1))

	p4 = plot_twotau_general(-evaps)
	plot!(p4,ylabel="evap + r_evap [g/kg]",ylims=(0,3),legend=false,yticks=collect(0:0.5:3),
		   left_margin=2mm,title="(d)")

	p5 = plot_twotau_general(wbf)
	plot!(p5,ylabel="wbf [g/kg]",ylims=(0,3),legend=false,yticks=collect(0:0.5:3),
		   title="(e)")
	
	p6 = plot_twotau_general(RH_tot)
	plot!(p6,ylabel="RH [%]",ylims=(75,120),legend=false,yticks=collect(80:10:120),
		   title="(f)")
	print("halway")
	
	p7 = plot_twotau_general(cnd)
	plot!(p7,ylabel="cond [g/kg]",ylims=(0,12),legend=false,yticks=collect(0:2:12),
		   title="(g)")
	
	p8 = plot_twotau_general(depo)
	plot!(p8,ylabel="depo [g/kg]",ylims=(0,12),legend=false,yticks=collect(0:2:12),
		   title="(h)")
	
	p9 = plot_twotau_general(qx_rim)
	plot!(p9,ylabel="qx_rim [g/kg]",ylims=(0,12),legend=false,yticks=collect(0:2:12),
		   title="(i)")
	
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=2)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=2,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=7)
	GC.gc()
	
	l9 = @layout [a{0.09h}
		       b c d
		       e f g
		       h i j]
	p_out = plot(blankp,p1,p2,p3,p4,p5,p6,p7,p8,p9,layout=l9,
		      dpi=:400,labelfontsize=6,tickfontsize=5,titlefontsize=6,
		      title_position=:left,size=(600,600),left_margin=5mm,bottom_margin=5mm)

	if save==true
                savefig(p_out,savepath*"micro_time.png")
        else
                return p_out
        end
end


function plot_hyds_time(save=false)
	qc = Array(dsnn["qc"]) .* 1000;
	qr = Array(dsnn["qr"]) .* 1000;
	qi = Array(dsnn["qi"]) .* 1000;
	qs = Array(dsnn["qs"]) .* 1000;
	qg = Array(dsnn["qg"]) .* 1000;
	qh = Array(dsnn["qh"]) .* 1000;
	
	qnc = Array(dsnn["qnc"]);
	qnr = Array(dsnn["qnr"]);
	qni = Array(dsnn["qni"]);
	qns = Array(dsnn["qns"]);
	qng = Array(dsnn["qng"]);
	qnh = Array(dsnn["qnh"]);

	p1 = plot_twotau_general(qc)
	plot!(p1,ylabel="qc [g/kg]",ylims=(0,0.33),legend=false,title="(a)")
	p2 = plot_twotau_general(qi)
	plot!(p2,ylabel="qi [g/kg]",ylims=(0,0.33),legend=false,title="(b)")
	p3 = plot_twotau_general(qs)
	plot!(p3,ylabel="qs [g/kg]",ylims=(0,0.33),legend=false,title="(c)")
	
	print("   one quarter    ")
	
	p4 = plot_twotau_general(qnc)
	plot!(p4,ylabel="qnc [1/kg]",legend=false,title="(d)")
	p5 = plot_twotau_general(qni)
	plot!(p5,ylabel="qni [1/kg]",legend=false,title="(e)")
	p6 = plot_twotau_general(qns)
	plot!(p6,ylabel="qns [1/kg]",legend=false,title="(f)")
	
	print("   two quarters    ")
	GC.gc()

	p7 = plot_twotau_general(qr)
	plot!(p7,ylabel="qr [g/kg]",ylims=(0,1.1),legend=false,title="(g)")
	p8 = plot_twotau_general(qg)
	plot!(p8,ylabel="qg [g/kg]",ylims=(0,1.1),yticks=collect(0:0.2:1.1),
	       legend=false,title="(h)")
	p9 = plot_twotau_general(qh)
	plot!(p9,ylabel="qh [g/kg]",ylims=(0,0.02),legend=false,title="(i)")
	
	print("   three quarters    ")
	GC.gc()
	
	p10 = plot_twotau_general(qnr)
	plot!(p10,ylabel="qnr [1/kg]",legend=false,title="(j)")
	p11 = plot_twotau_general(qng)
	plot!(p11,ylabel="qng [1/kg]",legend=false,title="(k)")
	p12 = plot_twotau_general(qnh)
	plot!(p12,ylabel="qnh [1/kg]",legend=false,title="(l)")
	
	blankp = plot(bg=:white,xaxis=false,yaxis=false,xticks=false,
		       yticks=false,xlim=(0,1),ylim=(0,1))
	
	plot!(blankp,[],c=:blue,fillrange=([NaN],[NaN]),fillalpha=0.4,label="convective",lw=0)
	plot!(blankp,[],c=:black,label="mean",lw=2)
	plot!(blankp,[],c=:coral,fillrange=([NaN],[NaN]),fillalpha=0.4,label="slow",lw=0)
	plot!(blankp,[],c=:black,label="median",lw=2,ls=:dash)
	plot!(blankp,legend=:top,legend_columns=2,legendfontsize=7)
	
	l12 = @layout [a{0.06h}
		       b c d
		       e f g
		       h i j
		       k l m]
	GC.gc()
	p_out = plot(blankp,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
		      layout=l12,dpi=:400,labelfontsize=7,tickfontsize=6,size=(700,800),
		      left_margin=2mm,titlefontsize=6,title_position=:left)
	if save==true
                savefig(p_out,savepath*"hyds_time.png")
        else
                return p_out
        end
end
