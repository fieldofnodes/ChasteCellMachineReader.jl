##################################################################
# Filename  : 20221222_add_general_data_to_csv_function.jl
# Author    : Jonathan Miller
# Date      : 2023-12-22
# Aim       : aim_script
#           : Add functionality to package to read any .dat file
#           : and turn to CSV
##################################################################
include("src/ChasteCellMachineReader.jl")
using .ChasteCellMachineReader

using DelimitedFiles
using DataFrames
using DataFramesMeta
using Chain


#https://github.com/fieldofnodes/ChasteDatFileToCSV.jl
sample_dir = "sample_data/20221222_encapsulated_sim_k1_1_single_realisation/"
data_names = ("cellscaling","cellages","machinedata","machinestate","cellorientation")
data_suffix = ".dat"
path = string.(sample_dir,data_names,data_suffix)

path = ("sample_data/20221222_encapsulated_sim_k1_1_single_realisation/cellscaling.dat",
    "sample_data/20221222_encapsulated_sim_k1_1_single_realisation/cellages.dat",
    "sample_data/20221222_encapsulated_sim_k1_1_single_realisation/machinedata.dat",
    "sample_data/20221222_encapsulated_sim_k1_1_single_realisation/machinestate.dat",
    "sample_data/20221222_encapsulated_sim_k1_1_single_realisation/cellorientation.dat")


data_df = get_cell_dataframe.(path)


