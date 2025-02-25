#---------------------------------------------------------------------------
#	"This script defines the functions that I use to conduct linear"
#	"interpolations of ICON data"
#       "Created by Cornelis Schwenk 06.10.2022"
#       "c.schwenk@uni-mainz.de"
#---------------------------------------------------------------------------

#	VITAL INFORMATION
#---------------------------------------------------------------------------
#	"The ICON data I am working with outputs data onto a lon/lat meshgrid"
#	"and onto model-levels, NOT pressure levels! These functions will"
#	"only work for this. Each coordinate has a fixed physical altitude"
#	"thoughout the simulation, which can be found in the icon_const file."
#	"If your ICON simulation writes onto pressure levels, this script will"
#	"NOT work for you!"
#---------------------------------------------------------------------------

#	NOTE
#---------------------------------------------------------------------------
#	"Some of the functions here are not very intuative. This is because"
#	"it was critical that the functions operate fast and efficiently and"
#	"with netCDF datasets, it is often not intuative what method of"
#	"extracting data is the fastest."
#---------------------------------------------------------------------------

#	INSTRUCTIONS
#---------------------------------------------------------------------------
#	"This script can only be included after the path to the icon data"
#	"AND the domain (1, 2 or 3 in my case) have been defined globally."
#	"Example:"
#			global DOMAIN = 2				 
#			global icon_path = "path/to/icon/output/directory/"
#			include("functions_lin_interp.jl")		 
#			interped_data = some_function(input1,input2,...)

#	"Once you have done this, you can only interpolate using the data"
#	"from that ICON simulation and the chosen domain. This is because"
#	"for efficiency reasons, global arrays such as LON, LAT and ALT are"
#	"defined."

#	"If you want to switch to a different ICON simulation or to a"
#	"different domain, you have to define the new DOMAIN and icon_path" 
#	"globally and include this script again, then execute functions"
#	"using Base.invokelatest."

#	"Example: you want to execute some_function(input1,input2,...) from"
#	"this script for DOMAIN3 and icon_run_1 and then for DOMAIN 1 and"
#	"icon_run_2. This is done as follows:"

#		global DOMAIN = 3
#		global icon_path = "path/to/output/icon/"
#		include("functions_lin_interp.jl")
#		interped_data_1 = Base.invokelatest(some_function,input1,input2,...)

#		global DOMAIN = 1
#		global icon_path = "path/to/output/icon/"
#		include("functions_lin_interp.jl")
#		interped_data_2 = Base.invokelatest(some_function,input1,input2,...)
		
#	"Overview of functions:"

#	"1) Interpolating a variable (e.g "temp") onto a lon, lat, pres and"  
#	"time point. [lon]=[lat]=°, [pres]=Pa, time as a DateTime object."
#	"Do not use this for multiple points using list comprehension. It"
#	"is VERY slow for that. See the next function instead."

#	interp_point_pres(lon,lat,pres,"temp",DateTime(2017,9,23,12,12))


#	"2) Interpolating a variable onto an array or lon, lat, pres and"
#	"time points. The time_array must be an array of DateTime points."
#	"This is very useful and fast for interpolating data onto a flight"
#	"path for instance. All arrays must ofcourse have same size."

#	interp_array_pres(lon_array,lat_array,pres_array,"variable",time_array)


#	"3) Interpolation onto an altitude (not pressure) point OR array."
#	"If you want to interpolate onto a single point, just use this also"
#	"and make the array with only one value like lon_array=[89.2] etc."

#	interp_array(lon_array,lat_array,z_array,"variable",time_array)	


#	"4) Interpolating a variable onto a pressure level meshgrid defined"
#	"by a lon_array and lat_array at a specific time. Here lon_array and"
#	"lat_array are 1dim arrays, plevel is a single Float and time a single"
#	"DateTime point. This can be used to upscale (or downscale) eulerian"
#	"output so that you can compare the output of simulations at different"
#	"scales."

#	interp_matrix_plevel(lon_array,lat_array,plevel,"variable",time_point)


#	"5) Interpolating onto a pressure level only, no interpolation in lon"
#	"and lat or time. So we are keeping the resolution as it is and just"
#	"finding what the output is on a pressure level. We are also not"
#	"interpolating in time, so time_point must be a time at which we wrote"
#	"ICON output."

#	"6) Writing the output of interp_matrix_plevel to a ncfile"
#	write_matrix(lon_array,lat_array,plevel,"temp",DateTime(2017,9,23,12))
#	"This creates a file called temp_plev20000_20170923_1200.nc"
#---------------------------------------------------------------------------

using Geodesy

#	BEGINNING OF SCRIPT
#---------------------------------------------------------------------------

#       create neccessary files
#--------------------------------
if (@isdefined DOMAIN) == false
	error("Error: DOMAIN needs to be defined first")
end
if (@isdefined icon_path) == false
	error("Error: icon_path needs to be defined first")
end

files = glob("mod_NWP*DOM0$(DOMAIN)*",icon_path)
const_file = glob("ICONCONST_DOM0$(DOMAIN)_DOM0$(DOMAIN)_ML_0001.nc",icon_path)
c_ds = Dataset(const_file)
global ALT = c_ds["z_mc"][:]
global dt_list = [DateTime(x[end-18:end-4],dateformat"yyyymmddTHHMMSS")+Millisecond(1) for x in files]
global LON = Dataset(files[1])["lon"][:]
global LAT = Dataset(files[1])["lat"][:]


#	define functions 
#--------------------------------
function nearest_t_neigh(dt)
	#"This function finds the two output files closest to our time"
        timear = abs.(dt_list .- dt)
        time_indx1 = findmin(timear)[2]
        timear[time_indx1] = timear[time_indx1]*100000000
        time_indx2 = findmin(timear)[2]

        p_t = sortperm([dt_list[time_indx1],dt_list[time_indx2]])
        return [time_indx1,time_indx2][p_t]
end

function nearest_neigh(lon,lat,z)
	#"This function finds the nearest lon, lat and altitude points"
	#"to our point. Imagine this as 8 points surrounding our desired"
	#"point like a box. Output is an array of LLA points. See Geodesy"
	#"for more info on LLA points https://github.com/JuliaGeo/Geodesy.jl"

        latar = abs.(LAT .- lat)
        lonar = abs.(LON .- lon)

        latindx1 = findmin(latar)[2]
        latar[latindx1] = 1000.0
        latindx2 = findmin(latar)[2]

        lonindx1 = findmin(lonar)[2]
        lonar[lonindx1] = 1000.0
        lonindx2 = findmin(lonar)[2]

        lat_inds = sort([latindx1,latindx2])
        lon_inds = sort([lonindx1,lonindx2])

        his = Vector{Int64}[]
        for lon_i in lon_inds
                for lat_i in lat_inds
                        ds_high = ALT[lon_i,lat_i,:]
                        highar = abs.(ds_high .- z)
                        hi1 = findmin(highar)[2]
                        #-------------------------------
                        # This makes sure the right heights are chosen
                        #-------------------------------
                        if z >= ds_high[hi1]
                                highar[hi1:end] .= 10000000
                        elseif z < ds_high[hi1]
                                highar[begin:hi1] .= 1000000
                        end
                        #--------------------------------
                        hi2 = findmin(highar)[2]

                        push!(his,sort([hi1,hi2]))
                end
        end

        points = [LLA(LAT[lat_inds[1]],LON[lon_inds[1]],ALT[lon_inds[1],lat_inds[1],his[1][1]]),
                   LLA(LAT[lat_inds[2]],LON[lon_inds[1]],ALT[lon_inds[1],lat_inds[2],his[2][1]]),
                   LLA(LAT[lat_inds[1]],LON[lon_inds[2]],ALT[lon_inds[2],lat_inds[1],his[3][1]]),
                   LLA(LAT[lat_inds[2]],LON[lon_inds[2]],ALT[lon_inds[2],lat_inds[2],his[4][1]]),
                   LLA(LAT[lat_inds[1]],LON[lon_inds[1]],ALT[lon_inds[1],lat_inds[1],his[1][2]]),
                   LLA(LAT[lat_inds[2]],LON[lon_inds[1]],ALT[lon_inds[1],lat_inds[2],his[2][2]]),
                   LLA(LAT[lat_inds[1]],LON[lon_inds[2]],ALT[lon_inds[2],lat_inds[1],his[3][2]]),
                   LLA(LAT[lat_inds[2]],LON[lon_inds[2]],ALT[lon_inds[2],lat_inds[2],his[4][2]])]

        return (points,(lon_inds[1],lon_inds[2],lat_inds[1],lat_inds[2]),his)
end

function ds_box(ds_var,lo1,lo2,la1,la2,his)
	#"This function returns a box of the values around the point."
	#"The logic is that you find the coordinates of the nearest neighbors"
	#"using the function above. You then pass the Dataset of your ICON"
	#"Data to this function and using the coordinates you obtain the box"
	#"of the values of the variable in question surrounding your point."
	ds_box1 = ds_var[lo1,la1,his[1][1]:his[1][2],1]
	ds_box3 = ds_var[lo1,la2,his[2][1]:his[2][2],1]
        ds_box5 = ds_var[lo2,la1,his[3][1]:his[3][2],1]
        ds_box7 = ds_var[lo2,la2,his[4][1]:his[4][2],1]
	return Float32[ds_box1[1],ds_box3[1],ds_box5[1],ds_box7[1],ds_box1[2],
                        ds_box3[2],ds_box5[2],ds_box7[2]]
end

function bilin_interp_alt(lo,la,z,indx,j)
	#"this function interpolates the ALTITUDE onto the plane above and/or"
	#"below our point. This sounds weird, but over mountainous terrain it"
	#"is possible that our model levels become very slanted. The four"
	#"corners of a data box are therefore at different altitudes and after"
	#"a bilinear interpolation of a variable it is neccessary to know onto"
	#"which altitude one has interpolated onto."

	points,(lo1,lo2,la1,la2),his = nearest_neigh(lo,la,z)
	vars = ds_box(ALT,lo1,lo2,la1,la2,his)

	y1 = -euclidean_distance(LLA(la,points[1].lon,z),LLA(points[1].lat,points[1].lon,z))
        x11 = -euclidean_distance(LLA(points[1].lat,points[1].lon,z),LLA(points[1].lat,lo,z))
        x12 = euclidean_distance(LLA(points[3].lat,points[3].lon,z),LLA(points[3].lat,lo,z))
	
        y2 = euclidean_distance(LLA(la,points[1].lon,z),LLA(points[2].lat,points[2].lon,z))
        x21 = -euclidean_distance(LLA(points[2].lat,points[2].lon,z),LLA(points[2].lat,lo,z))
        x22 = euclidean_distance(LLA(points[4].lat,points[4].lon,z),LLA(points[4].lat,lo,z))

        f_x_y1_z1 = ((x12 - 0)/(x12-x11))*vars[1] + ((0 - x11)/(x12-x11))*vars[3]
        f_x_y2_z1 = ((x22 - 0)/(x22-x21))*vars[2] + ((0 - x21)/(x22-x21))*vars[4]

        f_x_y_z1 = ((y2 - 0)/(y2-y1))*f_x_y1_z1 + ((0 - y1)/(y2-y1))*f_x_y2_z1

        f_x_y1_z2 = ((x12 - 0)/(x12-x11))*vars[5] + ((0 - x11)/(x12-x11))*vars[7]
        f_x_y2_z2 = ((x22 - 0)/(x22-x21))*vars[6] + ((0 - x21)/(x22-x21))*vars[8]

        f_x_y_z2 = ((y2 - 0)/(y2-y1))*f_x_y1_z2 + ((0 - y1)/(y2-y1))*f_x_y2_z2

	return (f_x_y_z1,f_x_y_z2)
end

function trilin_interp_lla_j(lo,la,z,vardict,dtime,indx,j)
	#"This function conducts a trilinear interpolation of the desired"
	#"variable onto our point. Think of vardict as a dataset."

	#"1) get the coordinates of the nearest neighbors and then the box."
	points,(lo1,lo2,la1,la2),his = nearest_neigh(lo,la,z)
	ds_var = vardict[indx[j]]
        vars = ds_box(ds_var,lo1,lo2,la1,la2,his)

	#"2) determine the altitude of the bilinear interpolated data on the"
	#"rectangular plane above and below our data point."
	z1,z2 = bilin_interp_alt(lo,la,z,indx,j)

	#"3) get all x and y distances from the point"
	y1 = -euclidean_distance(LLA(la,points[1].lon,z),LLA(points[1].lat,points[1].lon,z))
	x11 = -euclidean_distance(LLA(points[1].lat,points[1].lon,z),LLA(points[1].lat,lo,z))
	x12 = euclidean_distance(LLA(points[3].lat,points[3].lon,z),LLA(points[3].lat,lo,z))

	y2 = euclidean_distance(LLA(la,points[1].lon,z),LLA(points[2].lat,points[2].lon,z))
	x21 = -euclidean_distance(LLA(points[2].lat,points[2].lon,z),LLA(points[2].lat,lo,z))
	x22 = euclidean_distance(LLA(points[4].lat,points[4].lon,z),LLA(points[4].lat,lo,z))

	#"4) for the first plane do two linear interpolations in x"
	f_x_y1_z1 = ((x12 - 0)/(x12-x11))*vars[1] + ((0 - x11)/(x12-x11))*vars[3]
	f_x_y2_z1 = ((x22 - 0)/(x22-x21))*vars[2] + ((0 - x21)/(x22-x21))*vars[4]
	
	#"5) then interpolate between those points in y to complete the" 
	#"bilinear interpolation." 
	f_x_y_z1 = ((y2 - 0)/(y2-y1))*f_x_y1_z1 + ((0 - y1)/(y2-y1))*f_x_y2_z1

	#"6) do the same for the other plane"
	f_x_y1_z2 = ((x12 - 0)/(x12-x11))*vars[5] + ((0 - x11)/(x12-x11))*vars[7]
	f_x_y2_z2 = ((x22 - 0)/(x22-x21))*vars[6] + ((0 - x21)/(x22-x21))*vars[8]

	f_x_y_z2 = ((y2 - 0)/(y2-y1))*f_x_y1_z2 + ((0 - y1)/(y2-y1))*f_x_y2_z2

	#"7) interpolate in altitude between those two points to complete the"
	#"trilinear interpolation."
	f_xyz = f_x_y_z1 + (z-z1)*(f_x_y_z2-f_x_y_z1)/(z2-z1)
        return f_xyz
end

function interp_array(lo,la,z,var,dtime)
	#"This function finally interpolates onto a lon,lat,z and time point"
	#"or array by doing a final interpolation in time. For efficiency"
	#"reasons, the Datasets are loaded first into dictionaries and then"
	#"passed into the functions."
	
	if !allequal([length(lo),length(la),length(z),length(dtime)])
		error("All arrays lo, la, z and dtime must be of same length")
	end		

        k = 1
        counters = Int64[1]
        index1 = nearest_t_neigh(dtime[1])[1]
        indexn = nearest_t_neigh(dtime[end])[2]
        indices = [[x,x+1] for x in index1:(indexn-1)]

        vardict = Dict()
        for x in index1:indexn
                vardict[x] = Dataset(files[x])[var]
        end

        for dt in dtime
                if nearest_t_neigh(dt) != indices[k]
                        k+=1
                end
                push!(counters,k)
        end

        vals = Float64[]

        for i in 1:length(dtime)
                val1 = trilin_interp_lla_j(lo[i],la[i],z[i],vardict,dtime[i],indices[counters[i]],1)
                val2 = trilin_interp_lla_j(lo[i],la[i],z[i],vardict,dtime[i],indices[counters[i]],2)

                indx1 = indices[counters[i]][1]
                indx2 = indices[counters[i]][2]
                file1 = files[indx1]
                file2 = files[indx2]

                dt1 = DateTime(file1[end-18:end-4],dateformat"yyyymmddTHHMMSS")
                dt2 = DateTime(file2[end-18:end-4],dateformat"yyyymmddTHHMMSS")
                t1 = (dt1 - dtime[i]).value/1000
                t2 = (dt2 - dtime[i]).value/1000

                val = val1 + (0-t1)*(val2-val1)/(t2-t1)

                push!(vals,val)
        end
        return vals
end

function z_of_p(p1,p2,z1,z2,p)
	#"Im going to be honest, I don't think I use this function, but it"
	#"seems to determine the altitude of a pressure point between two points"
	#"directly above each other for which z1,z2 and p1,p2 are known."
	z = z1 .+ (log.(p) .- log(p1)) .* (z2-z1) ./ (log(p2) .- log(p1))
	return z
end

# "Now the same thing but with pressure, not altitude"
#--------------------------------------------------------------------------------
#"The functions below mostly do the same as above but instead using pressure"
#"coordinates in height rather than altitude. Pressure is a bit slower"
#"because we don't have a constant pressure at each ICON grid-box point"
#"throughout the simulation so I can't use a global P array like I do for ALT."

#"Since these are mostly the same functions, I have not commented them as"
#"explicitly."

function nearest_neigh_pres(lon,lat,p,pvardict,indx)

        latar = abs.(LAT .- lat)
        lonar = abs.(LON .- lon)

        latindx1 = findmin(latar)[2]
        latar[latindx1] = 1000.0
        latindx2 = findmin(latar)[2]

        lonindx1 = findmin(lonar)[2]
        lonar[lonindx1] = 1000.0
        lonindx2 = findmin(lonar)[2]

        lat_inds = sort([latindx1,latindx2])

        lon_inds = sort([lonindx1,lonindx2])
	
	P = pvardict[indx]

        ps = Vector{Int64}[]
        for lon_i in lon_inds
                for lat_i in lat_inds
                        ds_p = P[lon_i,lat_i,:,1]
                        p_ar = abs.(ds_p .- p)
                        p1 = findmin(p_ar)[2]
                        #-------------------------------
                        # This makes sure the right heights are chosen
                        #-------------------------------
                        if p >= ds_p[p1]
                                p_ar[p1:end] .= 10000000
                        elseif p < ds_p[p1]
                                p_ar[begin:p1] .= 10000000
                        end
                        #--------------------------------
                        p2 = findmin(p_ar)[2]

                        push!(ps,sort([p1,p2]))
                end
        end

        points = [LLA(LAT[lat_inds[1]],LON[lon_inds[1]],P[lon_inds[1],lat_inds[1],ps[1][1],1]),
                   LLA(LAT[lat_inds[2]],LON[lon_inds[1]],P[lon_inds[1],lat_inds[2],ps[2][1],1]),
                   LLA(LAT[lat_inds[1]],LON[lon_inds[2]],P[lon_inds[2],lat_inds[1],ps[3][1],1]),
                   LLA(LAT[lat_inds[2]],LON[lon_inds[2]],P[lon_inds[2],lat_inds[2],ps[4][1],1]),
                   LLA(LAT[lat_inds[1]],LON[lon_inds[1]],P[lon_inds[1],lat_inds[1],ps[1][2],1]),
                   LLA(LAT[lat_inds[2]],LON[lon_inds[1]],P[lon_inds[1],lat_inds[2],ps[2][2],1]),
                   LLA(LAT[lat_inds[1]],LON[lon_inds[2]],P[lon_inds[2],lat_inds[1],ps[3][2],1]),
                   LLA(LAT[lat_inds[2]],LON[lon_inds[2]],P[lon_inds[2],lat_inds[2],ps[4][2],1])]

        return (points,(lon_inds[1],lon_inds[2],lat_inds[1],lat_inds[2]),ps)
end
#------------------------------

function bilin_interp_lla_j_pres(lo,la,p,pvardict,indx,j,points,lo1,lo2,la1,la2,ps)
	#"This function interpolates the PRESSURE onto the plane above and below"
	#"This is because of the slant box problem."
	
	#points,(lo1,lo2,la1,la2),ps = nearest_neigh_pres(lo,la,p,pvardict,indx[j])
	
	P = pvardict[indx[j]]
	vars = ds_box(P,lo1,lo2,la1,la2,ps)

	y1 = -euclidean_distance(LLA(la,points[1].lon,log(p)),LLA(points[1].lat,points[1].lon,log(p)))
        x11 = -euclidean_distance(LLA(points[1].lat,points[1].lon,log(p)),LLA(points[1].lat,lo,log(p)))
        x12 = euclidean_distance(LLA(points[3].lat,points[3].lon,log(p)),LLA(points[3].lat,lo,log(p)))
	
        y2 = euclidean_distance(LLA(la,points[1].lon,log(p)),LLA(points[2].lat,points[2].lon,log(p)))
        x21 = -euclidean_distance(LLA(points[2].lat,points[2].lon,log(p)),LLA(points[2].lat,lo,log(p)))
        x22 = euclidean_distance(LLA(points[4].lat,points[4].lon,log(p)),LLA(points[4].lat,lo,log(p)))

        f_x_y1_z1 = ((x12 - 0)/(x12-x11))*vars[1] + ((0 - x11)/(x12-x11))*vars[3]
        f_x_y2_z1 = ((x22 - 0)/(x22-x21))*vars[2] + ((0 - x21)/(x22-x21))*vars[4]

        f_x_y_z1 = ((y2 - 0)/(y2-y1))*f_x_y1_z1 + ((0 - y1)/(y2-y1))*f_x_y2_z1

        f_x_y1_z2 = ((x12 - 0)/(x12-x11))*vars[5] + ((0 - x11)/(x12-x11))*vars[7]
        f_x_y2_z2 = ((x22 - 0)/(x22-x21))*vars[6] + ((0 - x21)/(x22-x21))*vars[8]

        f_x_y_z2 = ((y2 - 0)/(y2-y1))*f_x_y1_z2 + ((0 - y1)/(y2-y1))*f_x_y2_z2

	return (f_x_y_z1,f_x_y_z2)
end

function trilin_interp_lla_j_pres(lo,la,p,pvardict,vardict,dtime,indx,j)
	
	points,(lo1,lo2,la1,la2),ps = nearest_neigh_pres(lo,la,p,pvardict,indx[j])
	ds_var = vardict[indx[j]]
	vars = ds_box(ds_var,lo1,lo2,la1,la2,ps)
	
	p1,p2 = bilin_interp_lla_j_pres(lo,la,p,pvardict,indx,j,points,lo1,lo2,la1,la2,ps)

	y1 = -euclidean_distance(LLA(la,points[1].lon,log(p)),LLA(points[1].lat,points[1].lon,log(p)))
	x11 = -euclidean_distance(LLA(points[1].lat,points[1].lon,log(p)),LLA(points[1].lat,lo,log(p)))
	x12 = euclidean_distance(LLA(points[3].lat,points[3].lon,log(p)),LLA(points[3].lat,lo,log(p)))

	y2 = euclidean_distance(LLA(la,points[1].lon,log(p)),LLA(points[2].lat,points[2].lon,log(p)))
	x21 = -euclidean_distance(LLA(points[2].lat,points[2].lon,log(p)),LLA(points[2].lat,lo,log(p)))
	x22 = euclidean_distance(LLA(points[4].lat,points[4].lon,log(p)),LLA(points[4].lat,lo,log(p)))
	
	f_x_y1_z1 = ((x12 - 0)/(x12-x11))*vars[1] + ((0 - x11)/(x12-x11))*vars[3]
	f_x_y2_z1 = ((x22 - 0)/(x22-x21))*vars[2] + ((0 - x21)/(x22-x21))*vars[4]
	f_x_y_z1 = ((y2 - 0)/(y2-y1))*f_x_y1_z1 + ((0 - y1)/(y2-y1))*f_x_y2_z1

	f_x_y1_z2 = ((x12 - 0)/(x12-x11))*vars[5] + ((0 - x11)/(x12-x11))*vars[7]
	f_x_y2_z2 = ((x22 - 0)/(x22-x21))*vars[6] + ((0 - x21)/(x22-x21))*vars[8]
	f_x_y_z2 = ((y2 - 0)/(y2-y1))*f_x_y1_z2 + ((0 - y1)/(y2-y1))*f_x_y2_z2

	f_xyz = f_x_y_z1 + (log(p)-log(p1))*(f_x_y_z2-f_x_y_z1)/(log(p2)-log(p1))

        return f_xyz
end

function interp_array_pres(lo,la,p,var,dtime)
        k = 1
        counters = Int64[1]

        index1 = nearest_t_neigh(dtime[1])[1]
        indexn = nearest_t_neigh(dtime[end])[2]

        indices = [[x,x+1] for x in index1:(indexn-1)]

        vardict = Dict()
	pvardict = Dict()
        for x in index1:indexn
                vardict[x] = Dataset(files[x])[var]
		pvardict[x] = Dataset(files[x])["pres"][:]
        end

        for dt in dtime
                if nearest_t_neigh(dt) != indices[k]
                        k+=1
                end
                push!(counters,k)
        end

        vals = Float64[]

        for i in 1:length(dtime)
                val1 = trilin_interp_lla_j_pres(lo[i],la[i],p[i],pvardict,
						 vardict,dtime[i],indices[counters[i]],1)
                val2 = trilin_interp_lla_j_pres(lo[i],la[i],p[i],pvardict,
						 vardict,dtime[i],indices[counters[i]],2)

                indx1 = indices[counters[i]][1]
                indx2 = indices[counters[i]][2]
                file1 = files[indx1]
                file2 = files[indx2]

                dt1 = DateTime(file1[end-18:end-4],dateformat"yyyymmddTHHMMSS")
                dt2 = DateTime(file2[end-18:end-4],dateformat"yyyymmddTHHMMSS")
                t1 = (dt1 - dtime[i]).value/1000
                t2 = (dt2 - dtime[i]).value/1000

                val = val1 + (0-t1)*(val2-val1)/(t2-t1)

                push!(vals,val)
        end
	GC.gc()
        return vals
end

function interp_point_pres(lo,la,p,var,dtime)
	indx = nearest_t_neigh(dtime)
	vardict = Dict()
	pvardict = Dict()
	for x in indx[1]:indx[2]#
		vardict[x] = Dataset(files[x])[var]
		pvardict[x] = Dataset(files[x])["pres"]
	end

	val1 = trilin_interp_lla_j_pres(lo,la,p,pvardict,vardict,dtime,indx,1)
	val2 = trilin_interp_lla_j_pres(lo,la,p,pvardict,vardict,dtime,indx,2)

	indx1 = indx[1]
	indx2 = indx[2]
	file1 = files[indx1]
	file2 = files[indx2]

	dt1 = DateTime(file1[end-18:end-4],dateformat"yyyymmddTHHMMSS")
	dt2 = DateTime(file2[end-18:end-4],dateformat"yyyymmddTHHMMSS")
	t1 = (dt1 - dtime).value/1000
	t2 = (dt2 - dtime).value/1000

	val = val1 + (0-t1)*(val2-val1)/(t2-t1)

	return val
end

function interp_matrix_plevel(lo_ar::Vector{Float64},la_ar::Vector{Float64},p::Int,vr::String,dtime::DateTime)	
	if maximum(lo_ar) > maximum(LON)
		error("lon out of bounds, too large")
	elseif minimum(lo_ar) < minimum(LON)
		error("lon out of bounds, too small")
	elseif maximum(la_ar) > maximum(LAT)
		error("lat out of bounds, too large")
	elseif minimum(la_ar) < minimum(LAT)
		error("lat our of bounds, too small")
	end
		
	a = zeros(length(la_ar),length(lo_ar))
	
	for la in 1:length(la_ar)
		if la%10==0
			print("$(la)/$(length(la_ar)), ")
		end
		a[la,:] = interp_array_pres(lo_ar,[la_ar[la] for x in 1:length(lo_ar)],
					     [p for x in 1:length(lo_ar)],vr,
					     [dtime for x in 1:length(lo_ar)])
	end
	return Matrix{Float64}(a')
end


function write_matrix(lo_ar,la_ar,p,vr::String,dtime,dom::String)
	a = interp_matrix_plevel(lo_ar,la_ar,p,vr,dtime)
	dtstr = Dates.format(dtime, "yyyymmdd_HHMM")
	flnm = "ncfiles/$(vr)_plev$(p)_$(dtstr)_$(dom).nc"
	nccreate(flnm, vr, "lon", length(lo_ar), "lat", length(la_ar), 
		   Dict("units"=>"units"))
	ncwrite(a,flnm,vr)
end


#---------------------------------------------------------------------
#	"These functions don't interpolate in time"
#---------------------------------------------------------------------


function bilin_interp_lla_pres(lo,la,p,pvardict,indx,points,lo1,lo2,la1,la2,ps)
	P = pvardict[indx]
	vars = ds_box(P,lo1,lo2,la1,la2,ps)

	y1 = -euclidean_distance(LLA(la,points[1].lon,log(p)),LLA(points[1].lat,points[1].lon,log(p)))
        x11 = -euclidean_distance(LLA(points[1].lat,points[1].lon,log(p)),LLA(points[1].lat,lo,log(p)))
        x12 = euclidean_distance(LLA(points[3].lat,points[3].lon,log(p)),LLA(points[3].lat,lo,log(p)))
	
        y2 = euclidean_distance(LLA(la,points[1].lon,log(p)),LLA(points[2].lat,points[2].lon,log(p)))
        x21 = -euclidean_distance(LLA(points[2].lat,points[2].lon,log(p)),LLA(points[2].lat,lo,log(p)))
        x22 = euclidean_distance(LLA(points[4].lat,points[4].lon,log(p)),LLA(points[4].lat,lo,log(p)))

        f_x_y1_z1 = ((x12 - 0)/(x12-x11))*vars[1] + ((0 - x11)/(x12-x11))*vars[3]
        f_x_y2_z1 = ((x22 - 0)/(x22-x21))*vars[2] + ((0 - x21)/(x22-x21))*vars[4]

        f_x_y_z1 = ((y2 - 0)/(y2-y1))*f_x_y1_z1 + ((0 - y1)/(y2-y1))*f_x_y2_z1

        f_x_y1_z2 = ((x12 - 0)/(x12-x11))*vars[5] + ((0 - x11)/(x12-x11))*vars[7]
        f_x_y2_z2 = ((x22 - 0)/(x22-x21))*vars[6] + ((0 - x21)/(x22-x21))*vars[8]

        f_x_y_z2 = ((y2 - 0)/(y2-y1))*f_x_y1_z2 + ((0 - y1)/(y2-y1))*f_x_y2_z2

	return (f_x_y_z1,f_x_y_z2)
end

function trilin_interp_lla_pres(lo,la,p,pvardict,vardict,dtime,indx)

	points,(lo1,lo2,la1,la2),ps = nearest_neigh_pres(lo,la,p,pvardict,indx)
		
	ds_var = vardict[indx]

	vars = ds_box(ds_var,lo1,lo2,la1,la2,ps)

	p1,p2 = bilin_interp_lla_pres(lo,la,p,pvardict,indx,points,lo1,lo2,la1,la2,ps)

	y1 = -euclidean_distance(LLA(la,points[1].lon,log(p)),LLA(points[1].lat,points[1].lon,log(p)))
	x11 = -euclidean_distance(LLA(points[1].lat,points[1].lon,log(p)),LLA(points[1].lat,lo,log(p)))
	x12 = euclidean_distance(LLA(points[3].lat,points[3].lon,log(p)),LLA(points[3].lat,lo,log(p)))

	y2 = euclidean_distance(LLA(la,points[1].lon,log(p)),LLA(points[2].lat,points[2].lon,log(p)))
	x21 = -euclidean_distance(LLA(points[2].lat,points[2].lon,log(p)),LLA(points[2].lat,lo,log(p)))
	x22 = euclidean_distance(LLA(points[4].lat,points[4].lon,log(p)),LLA(points[4].lat,lo,log(p)))
	
	f_x_y1_z1 = ((x12 - 0)/(x12-x11))*vars[1] + ((0 - x11)/(x12-x11))*vars[3]
	f_x_y2_z1 = ((x22 - 0)/(x22-x21))*vars[2] + ((0 - x21)/(x22-x21))*vars[4]

	f_x_y_z1 = ((y2 - 0)/(y2-y1))*f_x_y1_z1 + ((0 - y1)/(y2-y1))*f_x_y2_z1

	f_x_y1_z2 = ((x12 - 0)/(x12-x11))*vars[5] + ((0 - x11)/(x12-x11))*vars[7]
	f_x_y2_z2 = ((x22 - 0)/(x22-x21))*vars[6] + ((0 - x21)/(x22-x21))*vars[8]

	f_x_y_z2 = ((y2 - 0)/(y2-y1))*f_x_y1_z2 + ((0 - y1)/(y2-y1))*f_x_y2_z2

	f_xyz = f_x_y_z1 + (log(p)-log(p1))*(f_x_y_z2-f_x_y_z1)/(log(p2)-log(p1))

        return f_xyz
end

function interp_array_pres_1time(lo,la,p,var,dtime)

	##	NO TIME INTERPOLATION HERE
	##	USE EXACT TIMES FROM ICON OUTPUT!!!

        index1 = argmin(abs.(dtime .- dt_list))

        vardict = Dict()
        pvardict = Dict()
               
	vardict[index1] = Dataset(files[index1])[var]
        pvardict[index1] = Dataset(files[index1])["pres"][:]

        vals = Float64[]

        for i in 1:length(lo)
                val = trilin_interp_lla_pres(lo[i],la[i],p[i],pvardict,
                                                 vardict,dtime,index1)

                push!(vals,val)
        end
        return vals
end


#	ONLY PRESSURE INTERPOLATION 
#	ONTO LON-LATS ALREADY THERE

function interp_var_domain_plev(plev,vr,dtime)
	index1 = argmin(abs.(dtime .- dt_list))
	ds = Dataset(files[index1])
	v = ds[vr][:]
	p = ds["pres"][:]
	
	var_mat = zeros(length(LON),length(LAT))
	
	for lo in 1:length(LON)
		for la in 1:length(LAT)
			p_col = p[lo,la,:]
			p_ar = abs.(plev .- p_col)
			
			pind1 = argmin(p_ar)
			p_ar[pind1] = 10000000
			pind2 = argmin(p_ar)
		
			pinds = sort([pind1,pind2])
			v1 = v[lo,la,pinds[1],1]
			v2 = v[lo,la,pinds[2],1]

			v_ = v1 + (log(plev)-log(p_col[pinds[1]]))*(v2-v1)/(log(p_col[pinds[2]])-log(p_col[pinds[1]]))
			var_mat[lo,la] = v_
		end
	end
	return var_mat
end

