#"This scripts writes a file that contains some important"
#"microphysical variables so that they don't have to be"
#"recalculated each time"

#"NOTE: you need the realtime file from zenodo and then"
#"      make the WCB_tau_vars_normed.jl file for this."

include("icon_WCB_functions.jl")

DR_time = calc_DR_time("nest")
coldf = sum([(Array(dsnn["q$(hy)_in"]) .- Array(dsnn["q$(hy)_out"])) for hy in ["i","s","h","g"]])
rainf = Array(dsnn["qr_in"]) .- Array(dsnn["qr_out"])
allfl = coldf .+ rainf
qhy = sum([Array(dsnn["q$(hy)"]) for hy in hymet])
DR_mix_time = calc_DR_mix_time("nest",qhy,allfl)
DR_mphys_time = calc_DR_mphys_time("nest",qhy,allfl)
PEt = calc_PE_time("nest",qhy,allfl)
CRt = calc_CR_time("nest")

ds_out = Dataset("../nest_data/PECR.nc","c")
defDim(ds_out,"time",101)
defDim(ds_out,"ntraj",length(t6n))

qhy_out = defVar(ds_out,"qhy",Float64,("ntraj","time"))
DR_out = defVar(ds_out,"DR",Float64,("ntraj","time"))
DRmix_out = defVar(ds_out,"DR_mix",Float64,("ntraj","time"))
DRmphys_out = defVar(ds_out,"DR_mphys",Float64,("ntraj","time"))
PE_out = defVar(ds_out,"PE",Float64,("ntraj","time"))
CR_out = defVar(ds_out,"CR",Float64,("ntraj","time"))
allfl_out = defVar(ds_out,"allfl",Float64,("ntraj","time"))
coldf_out = defVar(ds_out,"coldf",Float64,("ntraj","time"))
rainf_out = defVar(ds_out,"rainf",Float64,("ntraj","time"))

qhy_out[:,:] = qhy
DR_out[:,:] = DR_time
DRmix_out[:,:] = DR_mix_time
DRmphys_out[:,:] = DR_mphys_time
PE_out[:,:] = PEt
CR_out[:,:] = CRt
allfl_out[:,:] = allfl
coldf_out[:,:] = coldf
rainf_out[:,:] = rainf

close(ds_out)
