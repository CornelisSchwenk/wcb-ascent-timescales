---------------------------------------------------------------------------
#   "Skript for calculating the mean wind shear accumulated for all points"
#	"along the ascent as a function of tau_600. This script is sent to the"
#	"MOGON-II supercomputer by the do_calc_wind_shear.slurm file"

#       "Created by Cornelis Schwenk 22.12.2023"
#       "c.schwenk@uni-mainz.de"
#---------------------------------------------------------------------------

#	"This script uses the functions_lin_interp.jl, see that script for more"
#	"information. But suffice to say that you need to define the DOMAIN"
#	"and the icon path globally."

include("icon_WCB_functions.jl")
DOMAIN = 3
icon_path = "path/to/icon/ouput/"
include("functions_lin_interp.jl")
using Geodesy
using Base.Threads

#	"Calculate global variables that functions will use"
calc=true
if calc==true
	wcb,a,t = get_w_a_all();
	lor = Array(dsnr["lon"]) .* 180/pi;
	lar = Array(dsnr["lat"]) .* 180/pi;
	pr = Array(dsnr["p"]) ./ 100;
	zr = Array(dsnr["alt"]);
	t6n = Array(dsnn["tau_600"]) ./ 2;
	ds_dict = Dict()
	for r in 0.5:0.5:96
		ds_dict[r] = Dataset(nest_NWP_fls_dom3[Int(r*4) + 1])
	end
	for i in 1:400310
		if a[i] > 1
			lor[i,1:a[i]] .= NaN
			lar[i,1:a[i]] .= NaN
			pr[i,1:a[i]] .= NaN
			zr[i,1:a[i]] .= NaN
			
			lor[i,a[i] + Int(t6n[1] .* 2):end] .= NaN
			lar[i,a[i] + Int(t6n[1] .* 2):end] .= NaN
			pr[i,a[i] + Int(t6n[1] .* 2):end] .= NaN
			zr[i,a[i] + Int(t6n[1] .* 2):end] .= NaN
		end
	end
end

#	"Define what min and maximum Plevel you want to look at"
PLEV = [600,1030];

#	"Get indices of nearest points"
function get_nearest(rt,ix::Int)
	lo = lor[ix,Int(rt*2)]
	la = lar[ix,Int(rt*2)]
	z = zr[ix,Int(rt*2)]
	return nearest_neigh(lo,la,z)
end

function get_wind_shear_Z(rt,ix::Int,dsu,dsv)
	points,(lo1,lo2,la1,la2),his = get_nearest(rt,ix)
	u = ds_box(dsu,lo1,lo2,la1,la2,his)
	v = ds_box(dsv,lo1,lo2,la1,la2,his)
	z = ds_box(ALT,lo1,lo2,la1,la2,his)

	du = u[1:4] .- u[5:end]
	dv = v[1:4] .- v[5:end]
	dz = z[1:4] .- z[5:end]
	uv = sqrt.((du .^ 2) .+ (dv .^ 2))
	wind_shear = mean(uv ./ dz)
	return wind_shear
end

function get_wind_shear_SN(rt,ix::Int,dsu,dsv)
	points,(lo1,lo2,la1,la2),his = get_nearest(rt,ix)
	u = ds_box(dsu,lo1,lo2,la1,la2,his)
	v = ds_box(dsv,lo1,lo2,la1,la2,his)
	z = ds_box(ALT,lo1,lo2,la1,la2,his)

	du = u[1:2:end] .- u[2:2:end]
	dv = v[1:2:end] .- v[2:2:end]
	point1 = LLA(LAT[la1],LON[lo1],mean(z))
	point2 = LLA(LAT[la2],LON[lo1],mean(z))
	dx = euclidean_distance(point1,point2) / 1000 

	uv = sqrt.((du .^ 2) .+ (dv .^ 2))
	wind_shear = mean(uv ./ dx)
	return wind_shear
end

function get_wind_shear_EW(rt,ix::Int,dsu,dsv)
	points,(lo1,lo2,la1,la2),his = get_nearest(rt,ix)
	u = ds_box(dsu,lo1,lo2,la1,la2,his)
	v = ds_box(dsv,lo1,lo2,la1,la2,his)
	z = ds_box(ALT,lo1,lo2,la1,la2,his)
	du = [(u[1] - u[3]),(u[2] - u[4]),(u[5] - u[7]),(u[6] - u[8])]
	dv = [(v[1] - v[3]),(u[2] - v[4]),(v[5] - v[7]),(v[6] - v[8])]
	point1 = LLA(LAT[la1],LON[lo1],mean(z))
	point2 = LLA(LAT[la1],LON[lo2],mean(z))
	dx = euclidean_distance(point1,point2) / 1000 

	uv = sqrt.((du .^ 2) .+ (dv .^ 2))
	wind_shear = mean(uv ./ dx)
	return wind_shear
end


function get_all_shear_t6(plev,t6,func)
	wz_ar = []
	@threads for r in 0.5:0.5:96
		inds = findall((t6n .== t6) .& (minlon3 .< lor[:,Int(r * 2)] .< maxlon3) .& 
					   (minlat3 .< lar[:,Int(r * 2)] .< maxlat3) .& 
					   (plev[1] .< pr[:,Int(r * 2)] .< plev[2]))
		if length(inds) .> 0
			dsu = ds_dict[r]["u"][:,:,:,1]
			dsv = ds_dict[r]["v"][:,:,:,1]
			push!(wz_ar,[func(r,i,dsu,dsv) for i in inds]...)
		
		end
		if r in [10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0] 
			GC.gc()
		end
	end
	return wz_ar
end

#----------WIND SHEAR-------------
println("mean EW wind shear for t = 1:0.5:48")
print([mean(get_all_shear_t6(PLEV,t,get_wind_shear_EW)) for t in 1:0.5:48])
println("----------------------------")
println("mean SN wind shear for t = 1:0.5:48")
print([mean(get_all_shear_t6(PLEV,t,get_wind_shear_SN)) for t in 1:0.5:48])
println("----------------------------")
println("mean Z wind shear for t = 1:0.5:48")
print([mean(get_all_shear_t6(PLEV,t,get_wind_shear_Z)) for t in 1:0.5:48])

