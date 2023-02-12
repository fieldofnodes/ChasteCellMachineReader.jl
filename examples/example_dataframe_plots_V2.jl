"""
    Sample script to take a folder of chaste data
    create a data frame of the time series for each simulation
    create some plots
    save to file
"""

using Pkg
Pkg.activate(".")
using ChasteCellMachineReader
using DataFrames
using DataFramesMeta
using Find
using CSV



t6ss_data_root = "/Users/jmille15/Data/chaste/T6SS/"
sim_folder = t6ss_data_root*"20221220_K1-10_Realisations-50"
find_pattern = "machinestate.dat"



"""
    Get cell paths
"""
cell_path = @chain sim_folder begin
    find(_,0,7,find_pattern) 
    filter(x->occursin("plots",x),_)
end

"""
    Choose which cell simulation to execute
"""
root_path = @chain cell_path begin
   dirname.(_) 
end

core_file_name = @chain cell_path begin
    dirname.(_) 
    dirname.(_) 
    basename.(_) 
 end

 



"""
    Generate and write to file each TS per chaste file
"""
output_desc = "TS_df"
file_ext = "csv"

output_dfs = string.(root_path,"/",
    output_desc,"_",
    core_file_name,".",file_ext)

cdfs = get_cell_dataframe_TS.(cell_path)
CSV.write.(output_dfs,cdfs)



"""
    Save plots to file
"""
plots = plot_cell_populations_TS.(cdfs)

each_plot = "plot_".*[
    "cell_populations",
    "machines_diff_state",
    "per_cap_machines_state",
    "total_machines",
    "per_cap_total_machines",
    "total_neighbours",
    "per_cap_total_neighbours"]

file_ext = "png"

output_dfs = [string.(root_path[r],"/",
    each_plot,"_",core_file_name[r],".",file_ext) 
    for r in eachindex(root_path)]


[save.(output_dfs[i],plots[i]) for i in eachindex(plots)]


