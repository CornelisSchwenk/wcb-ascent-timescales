using NaNMath; nm=NaNMath
using StatsBase, LsqFit, Glob, NCDatasets, Plots

#---------------------------------------------------------------------------
#	"Skript for the computation of thermodynamic variables etc."
#	"Created by Cornelis Schwenk 19.10.2023"
#	"c.schwenk@uni-mainz.de"
#---------------------------------------------------------------------------

global M_air = 28.96470 #"g/mol"
global M_H2O = 18.01528 #"g/mol"

function convert_specific_humidity(qv_ppm,T,rho) 
	#"converts ppm to g/kg"
	
	qv_ppm = prepare_input(qv_ppm)
	if qv_ppm < 0
		qv_ppm = NaN
	end
	T = prepare_input(T)
	rho = prepare_input(rho)
	
	# "Constants"
	NA = 6.02214076e23  # "Avogadro's number in molecules/mol"
	MW_water = 18.01528  # "Molar mass of water in g/mol"
	MW_air = 28.9647  # "Molar mass of air in g/mol"

	# "Convert qv from ppm to mole fraction"
	qv_mole_fraction = qv_ppm * 1e-6
	X = qv_mole_fraction
	specific_humidity = qv_x(X)
	return specific_humidity
end


function calculate_pressure_from_rho_and_T(rho, T)
	#"rho must be in Molecules per cubic meter"
	#"T must be in Kelvin"
	#"Boltzmann constant in J/K"
	k_B = 1.380649e-23

	# "Calculate pressure using the ideal gas law"
	P = rho .* k_B .* T
	return P
end

function qv_x(X) 
	#-----------------------------
	#"This function calculates the specific humidity in kg/kg"
	#"given the mole fraction of water vapor in air: X in [0,1]"
	#-----------------------------
	return X .* M_H2O ./ ((X .* M_H2O) .+ (1 .- X) .* M_air)
end

function qv_ppmv(p)
	#-----------------------------
	#"This function gives the specific humidity in kg/kg"
	#"given the volume fraction in ppmv"
	#-----------------------------
	return p .* M_H2O ./ (1e6 .* M_air)
end

function e_sat_i(T)
        #------------------------------
        #"This function calculates the saturation vapor pressure over ICE"
	#"T must be in Kelvin, result is given in Pa"
        #------------------------------
        if ismissing(T)
		T = NaN
	end
        if (T > 273.15) | (T < 110)
                return NaN
        else
                return exp(9.550426 - 5723.265/T + 3.53068*log(T) - 0.00728332*T)
	end
end

function e_sat_w(T)
        #------------------------------
        #"This function calculates the saturation vapor pressure over WATER"
	#"T must be in Kelvin, result is given in Pa"
        #------------------------------
        if ismissing(T)
		T = NaN
	end
	if (T > 332) | (T < 123)
                error("T out of bounds. T in (123,332)K")
        end
        return exp(54.842763 - (6763.22 / T) - 4.210*log(T) + 0.000367*T +
                   tanh(0.0415 * (T - 218.8))*(53.878 - (1331.22 / T) -
                                               9.44523*log(T) + 0.014025*T))
end

function qv_sat_w(T,p)
        #------------------------------
        #"This function calculates the saturation qv(T,p) over WATER"
	#"T must be in Kelvin, p in Pa and result is given in kg/kg"
        #------------------------------
        return e_sat_w(T)*0.622/p
end

function qv_sat_i(T,p)
        #------------------------------
        #"This function calculates the saturation qv(T,p) over ICE"
	#"T must be in Kelvin, p in Pa and result is given in kg/kg"
        #------------------------------
        return e_sat_i(T)*0.622/p
end

function RH_w(T,qv,p)
        #------------------------------
        #"This function calculates RH over WATER"
	#"T must be in Kelvin, p in Pa and qv in kg/kg"
        #"Result is given in %"
	#------------------------------
        e = (qv .* p) ./ 0.622
        return 100*e ./ e_sat_w.(T)
end

function RH_i(T,qv,p)
        #------------------------------
        #"This function calculates RH over ICE"
	#"T must be in Kelvin, p in Pa and qv in kg/kg"
        #"Result is given in %"
	#------------------------------
        e = (qv .* p) ./ 0.622
        return 100*e ./ e_sat_i.(T)
end

function Theta(T,p)
        #------------------------------
        #"This function calculates Potential Temperature"
	#"T must be in Kelvin, p in Pa"
        #"Result is given in Kelvin"
	#------------------------------
	return T .* ((100000 ./ p) .^ 0.286)
end
