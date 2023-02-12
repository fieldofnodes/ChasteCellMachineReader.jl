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
radius_k1_100_data = t6ss_data_root*"K1-100_CircularDomainSurroundedByAttackers"
range_k1_data = t6ss_data_root*"RangeOverK1_CircularDomainSurroundedByAttackers"
find_pattern = "machinestate.dat"


new_range_radius_k1_100_folder = t6ss_data_root*"RangeOverRadius_K1_100"
new_range_k1_folder = t6ss_data_root*"RangeOverK1"

"""
    Make directory for output files
"""
mkdir(new_range_radius_k1_100_folder)
mkdir(new_range_k1_folder)


"""
    Make directory for     
    Range over radii when K1 = 100
""" 
@chain radius_k1_100_data begin
    find(_,0,5,find_pattern) 
    new_folder_name.(_) 
    new_range_radius_k1_100_folder .* "/" .* _ 
    mkdir.(_)
end

"""
    Make directory for 
    Range over K1
""" 
@chain range_k1_data begin
    find(_,0,5,find_pattern) 
    new_folder_name.(_) 
    new_range_k1_folder .* "/" .* _ 
    mkdir.(_)
end



"""
    Get cell paths
"""
cell_path = @chain range_k1_data begin
    find(_,0,5,find_pattern) 
end

"""
    Choose which cell simulation to execute
"""
root_path = new_range_k1_folder





"""
    Copy chaste output cell file to new folder locations
"""
file_ext = "dat"
output_desc = "copy"

output_dfs = make_output_path.(root_path,cell_path, 
    output_desc,file_ext)

cp.(cell_path,output_dfs)
    

"""
    Generate and write to file each TS per chaste file
"""
output_desc = "TS_df"
file_ext = "csv"
output_dfs = make_output_path.(root_path,cell_path, 
    output_desc,file_ext)
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

output_dfs = [make_output_path.(
    root_path,c,each_plot,file_ext) 
    for c in cell_path]


[save.(output_dfs[i],plots[i]) for i in eachindex(plots)]