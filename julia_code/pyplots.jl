using PyPlot, PyCall, NCDatasets, Glob, Dates
using Images, ImageView

include("use_config.jl")
plt = pyimport("matplotlib.pyplot")
ccrs = pyimport("cartopy.crs")
mpl = pyimport("matplotlib")

#--------------------------------------
#	"define paths"
#--------------------------------------

nest_NWP_path = "/path/to/icon/output/"

Era_PATH = "/path/to/era5/data/"

savepath = "../plots/"
#--------------------------------------
#	"define necessary functions"
#--------------------------------------

function mesh(lon_,lat_)
        return [ones(length(lat_),length(lon_)) .* lon_',ones(length(lat_),length(lon_)) .* lat_]
end

function get_dsNWP(dt::DateTime)
        dtstr = Dates.format(dt,"yyyymmddTHHMMSS")
        dsNWP = Dataset(nest_NWP_path * "mod_NWP_global_DOM01_" * dtstr * "Z.nc")
        return dsNWP
end

function get_NWPvar(dt::DateTime,vr::String)
        dtstr = Dates.format(dt,"yyyymmddTHHMMSS")
        dsNWP = Dataset(nest_NWP_path * "mod_NWP_global_DOM01_" * dtstr * "Z.nc")
        out = Array(dsNWP[vr])
	close(dsNWP)
	return out
end

function get_Bdsera(dt)
        dtstr = Dates.format(dt,"yyyymmdd_HH")
        return Dataset(Era_PATH * "B$(dtstr).nc")
end

function get_Bdsera_var(dt,vr)
        dtstr = Dates.format(dt,"yyyymmdd_HH")
        dsera = Dataset(Era_PATH * "B$(dtstr).nc")
        return reverse(dsera[vr][451:831,21:261,1],dims=2)
        close(dsera)
end

#--------------------------------------
#	"load and create necessary data"
#--------------------------------------

ds_const_dom1 = Dataset(nest_NWP_path * "ICONCONST_DOM01_DOM01_ML_0001.nc")
lo1 = Array(ds_const_dom1["lon"])
la1 = Array(ds_const_dom1["lat"])
lo2d, la2d = mesh(lo1, la1)

ds_const_dom2 = Dataset(nest_NWP_path * "ICONCONST_DOM02_DOM02_ML_0001.nc")
lo2 = Array(ds_const_dom2["lon"])
la2 = Array(ds_const_dom2["lat"])

ds_const_dom3 = Dataset(nest_NWP_path * "ICONCONST_DOM03_DOM03_ML_0001.nc")
lo3 = Array(ds_const_dom3["lon"])
la3 = Array(ds_const_dom3["lat"])

#"for era5 data we only look at our DOMAIN1"
dsera1 = get_Bdsera(DateTime(2017,9,20,0))
eralo = dsera1["longitude"][451:831]
erala = reverse(dsera1["latitude"][21:261])
lo2dera,la2dera = mesh(eralo,erala)

#--------------------------------------
#	"define plotting functions"
#--------------------------------------

function case_study_plot(dt::DateTime,glob_nest::String)
        #"------load data--------"
        pmsl = get_NWPvar(dt,"pres_msl",glob_nest)[:,:,1] ./ 100
        temp = get_NWPvar(dt,"temp",glob_nest)[:,:,end,1] .- 273.15
        u = get_NWPvar(dt,"u",glob_nest)[:,:,end,1]
        v = get_NWPvar(dt,"v",glob_nest)[:,:,end,1]

        vmax,vmin = 25,0
        #"-----create figure and customize-----"
        fig = plt.figure(figsize=(10,8))

        ax1 = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
        ax1.set_extent([-65,5,30,80],
                       crs=ccrs.PlateCarree())
        ax1.coastlines(resolution="50m",color="white",lw=1.5)

        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                        linewidth=0.5,color="black", alpha=1, linestyle="--")

        bbox_ax = ax1.get_position()
        cb_ax_ver = fig.add_axes([0.95,bbox_ax.y0+0.01, 0.015, bbox_ax.y1-bbox_ax.y0 - 0.03])

        #"-------heatmap of temperature-------"
        qcs1 = ax1.contourf(lo2d,la2d, temp',
                        levels = LinRange(minimum(temp),maximum(temp),100),
                        cmap = "plasma",vmax=vmax,vmin=vmin)

        #"-------contourplot of sea-level pressure-------"
        qcs2 = ax1.contour(lo2d,la2d, pmsl',
                        levels = 800:5:1040,
                        colors="cyan",linewidths=1)

        ax1.clabel(qcs2, inline=1, fontsize=10)

        #"-------quiver plot of wind-------"
        skp=6
        quiv = ax1.quiver(lo2d[1:skp:end,1:skp:end],la2d[1:skp:end,1:skp:end],
                           u[1:skp:end,1:skp:end]',v[1:skp:end,1:skp:end]',
                    scale_units="xy", scale=9, width=0.0011)

        qk = plt.quiverkey(quiv,1,1.02,15,"15m/s",labelpos="E",coordinates="axes")

        #"-----colorbar------"
        sm1 = plt.cm.ScalarMappable(cmap="plasma",norm=qcs1.norm)
        sm1._A = []

        cbar1 = fig.colorbar(sm1,cax=cb_ax_ver,orientation="vertical",extend="both")
        cbar1.set_label("Temperature [°C]",size="xx-large")

        #"-----title-----"
        ax1.set_title("a) \n",fontsize = 18, loc="left")
        ax1.set_title("\n $(dt)",fontsize = 15, loc="center")

        dtstr = Dates.format(dt,"yyyymmddTHHMM")
        fig.savefig("$(savepath)$(dtstr)_$(glob_nest).png",dpi = 250,bbox_inches="tight")
        plt.close()
end

#---------------------------------------------
#	"Pressure level data plotting"
#---------------------------------------------

#"This uses the linear interpolation functions from functions_lin_interp.jl"

function case_study_PV_plev_plot(dt,plev,vmin,vmax,glob_nest::String) #"plev in hPa, vmin,vmax in PVU"
        #"------load data--------"

        global DOMAIN = 1
        global icon_path = pathdict[glob_nest]
        include("functions_lin_interp.jl")

        pv = Base.invokelatest(interp_var_domain_plev,plev*100,
                                "pv",dt) .* 1e6
        u = Base.invokelatest(interp_var_domain_plev,plev*100,
                                "u",dt)
        v = Base.invokelatest(interp_var_domain_plev,plev*100,
                                "v",dt)

        pmsl = get_NWPvar(dt,"pres_msl",glob_nest)[:,:,1] ./ 100

        #vmax,vmin = 9,0
        #"-----create figure and customize-----"
        fig = plt.figure(figsize=(10,8))

        ax1 = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
        ax1.set_extent([-65,5,30,80],
                       crs=ccrs.PlateCarree())
        ax1.coastlines(resolution="50m",color="k",lw=1)

        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                        linewidth=0.5,color="black", alpha=1, linestyle="--")

        bbox_ax = ax1.get_position()
        cb_ax_ver = fig.add_axes([0.95, bbox_ax.y0+0.01, 0.015, bbox_ax.y1-bbox_ax.y0 - 0.03])

        #"-------heatmap of temperature-------"
        qcs1 = ax1.contourf(lo2d,la2d, pv',
                        levels = LinRange(minimum(pv),maximum(pv),100),
                        cmap = "Wistia_r",vmax=vmax,vmin=vmin)

        #"-------contourplot of sea-level pressure-------"
        qcs2 = ax1.contour(lo2d,la2d, pmsl',
                        levels = 800:5:1040,
                        colors="blue",linewidths=1)

        ax1.clabel(qcs2, inline=1, fontsize=12)
        #"-------quiver plot of wind-------"
        skp=8
        quiv = ax1.quiver(lo2d[1:skp:end,1:skp:end],la2d[1:skp:end,1:skp:end],
                           u[1:skp:end,1:skp:end]',v[1:skp:end,1:skp:end]',
                    scale_units="xy", scale=20, width=0.0011)

        qk = plt.quiverkey(quiv,1,1.02,35,"35m/s",labelpos="E",coordinates="axes")

        #"-----colorbar------"
        sm1 = plt.cm.ScalarMappable(cmap="Wistia_r",norm=qcs1.norm)
        sm1._A = []

        cbar1 = fig.colorbar(sm1,cax=cb_ax_ver,orientation="vertical",extend="both")
        cbar1.set_label("Potential Vorticity [PVU]",size="xx-large")

        #"-----title-----"
        ax1.set_title("b) \n",fontsize = 18, loc="left")
        ax1.set_title("\n $(dt) @ $(plev) hPa",fontsize = 15, loc="center")
        GC.gc()
        dtstr = Dates.format(dt,"yyyymmddTHHMM")
        fig.savefig("$(savepath)PV_$(dtstr)_$(glob_nest).png",dpi = 250,bbox_inches="tight")
        plt.close()
end

#----------------------------------------------
#	"combining plots"
#----------------------------------------------

function comb_plot(dt,filename)
        dtstr = Dates.format(dt,"yyyymmddTHHMM")
        img1 = load(savepath * dtstr*"_glob.png")
        img2 = load(savepath * "PV_" * dtstr * "_glob.png")
        comb_img = hcat(img1, img2)
        save(savepath * "$(filename).png",comb_img)
end

function make_all_case_study_plots()
        dt_list = [DateTime(2017,9,20,8),DateTime(2017,9,21,4),DateTime(2017,9,22,4),
                    DateTime(2017,9,23,8),DateTime(2017,9,24,0)]
        nm_list = ["A1","A2","A3","A4","A5"]

        for dt in dt_list
                println(dt)
                case_study_plot(dt,"glob")
                case_study_PV_plev_plot(dt,200,0,9,"glob")
        end
        println("combining")
        for (i,dt) in enumerate(dt_list)
                comb_plot(dt,nm_list[i])
        end
end

#----------------------------------------------
#	"era5 data plotting"
#----------------------------------------------

function plot_era5_map(dt)
        dsera = get_Bdsera(dt)
	pmsl = get_Bdsera_var(dt,"msl") ./ 100
	tcc = get_Bdsera_var(dt,"tcc")

        fig = plt.figure(figsize = (10,8))

        ax1 = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
        ax1.set_extent([minimum(eralo),maximum(eralo),minimum(erala),maximum(erala)],
		       crs=ccrs.PlateCarree())
        ax1.coastlines(resolution="50m",color="k")
        
	gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=true, 
			linewidth=0.5,color="black", alpha=1, linestyle="--")

	bbox_ax = ax1.get_position()

	cb_ax_ver1 = fig.add_axes([0.95,bbox_ax.y0, 0.015, bbox_ax.y1-bbox_ax.y0])
	cb_ax_ver2 = fig.add_axes([1.03,bbox_ax.y0, 0.015, bbox_ax.y1-bbox_ax.y0])
	
	qcs1 = ax1.contourf(lo2dera,la2dera, tcc',
                        levels = LinRange(0,1,20), cmap = "Blues_r")
	
	qcs2 = ax1.contour(lo2dera,la2dera, pmsl',
			levels = LinRange(960,1030,20),cmap = "plasma")

	ax1.plot([minimum(lo2),maximum(lo2),maximum(lo2),minimum(lo2),minimum(lo2)],
		 [minimum(la2),minimum(la2),maximum(la2),maximum(la2),minimum(la2)],
		 c="red",linewidth=1.5)
	ax1.plot([minimum(lo3),maximum(lo3),maximum(lo3),minimum(lo3),minimum(lo3)],
		 [minimum(la3),minimum(la3),maximum(la3),maximum(la3),minimum(la3)],
		 c="black",linewidth=1.5)
	ax1.plot([-60,-20,-20,-60,-60],[30,30,50,50,30],c="black",linestyle="dashed",
		 linewidth=2)
	
	sm1 = plt.cm.ScalarMappable(cmap="Blues_r",norm=qcs1.norm)
	sm1._A = []
	
	sm2 = plt.cm.ScalarMappable(cmap="plasma",norm=qcs2.norm)
	sm2._A = []
	
	cbar1 = fig.colorbar(sm1,cax=cb_ax_ver1,orientation="vertical")
	cbar1.set_label("Total Cloud Cover []",size="x-large")	
	
	cbar2 = fig.colorbar(sm2,cax=cb_ax_ver2,orientation="vertical",extend="both")
	cbar2.set_label("Sea Level Pressure [hPa]",size="x-large")	
	
	#return fig
        dtstr = Dates.format(dt,"yyyymmddTHHMM")
	svpth = "../plots/case_study/"
        fig.savefig("$(savepath)Fig1.png",dpi = 250,bbox_inches="tight")
	plt.close()
end
	
#-----------------------------------------------
#	"Execute plotting functions"
#-----------------------------------------------
function main()
	make_all_case_study_plots
	plot_era5_map(DateTime(2017, 9, 23, 10), save_path)
end

#"This will produce the plots and they will be chronological if sorted by date."
#"This means that the first one will be FigA1, the second FigA2 and so on."
main()

