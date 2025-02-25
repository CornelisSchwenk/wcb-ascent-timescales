#---------------------------------------------------------
#	"Created by Cornelis Schwenk 12.12.2023 c.schwenk@uni-mainz.de"
#
#	"!! You can only do this AFTER you have download the files from"
#	"10.5281/zenodo.12705854 and put them into nest_data/ !!"
#---------------------------------------------------------

include("WCB_selection.jl")

#"print date and time just for completeness"

dt_start = Dates.now()
println(Dates.format(dt_start, "HH:MM:SS"))

write_tauWCB_var_file_normed()

#"Print time and date. You can do this if you want to know how long it all took"
dt_end = Dates.now()
println(Dates.format(dt_end, "HH:MM:SS"))
println("done")
println("this took $(round(dt_end - dt_start,Minute))")
