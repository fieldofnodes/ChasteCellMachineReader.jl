module ChasteCellMachineReader

using Reexport

@reexport using Revise
@reexport using Chain
@reexport using DataFrames
@reexport using DataFramesMeta
@reexport using DelimitedFiles
@reexport using CairoMakie
@reexport using Dates
@reexport using Statistics


export
    MachineState,
    MachineData,
    MachineStateCellProperties,
    MachineDataProperties,
    MachineStateTimeSeries,
    PopulationCountsPerTimeStep,
    TimeVecVariants,
    ApoptosisTime,
    Radius,
    Cylinder,
    Capsule,
    N₀AlreadySet,
    CellLabelsTimeSeries,
    CellTimeSeries,
    TimeForDeathRate,
    DeathRate,
    DeathRateTime,
    get_cell_data,
    get_cell_dataframe,
    shift_time_back_to_zero,
    get_time_vec_variants,
    get_type_count_per_time,
    get_population_counts_per_time_step,
    get_competition_score_function,
    get_population_competition_time_series,
    cₐ,
    cₗ,
    today_nog_gap,
    N₀,
    R₀,
    dᵣ,
    β,
    get_norm_vec,
    shift_scale_uniform,
    get_mean_cycle_time,
    get_time_series_labels_diff_to_cell,
    get_cell_population_ts,
    get_count_neighbour_targets,
    get_mean_neighbours_at_time_time_slices,
    get_mean_competition,
    compute_target_death_rate,
    get_competition_per_params,
    get_mean_comp,
    plot_population_neighbours_competition,
    convert_mat_to_df,
    compute_lm_predict,
    fit_line,
    compute_pop_ode_sol,
    compute_time_extinction,
    simulation_population_TS,
    simulation_analytical,
    get_neighbour_radius,
    get_radius,
    get_realisation_number,
    get_Nₜ,
    get_Nₐ,
    get_Nt,
    add_realisations,
    get_mean_pop_with_std_errors,
    std_error,
    missing_to_zero,
    convert_missing_to_zero,
    if_missing_then_zero,
    input_neg_lookbehind,
    parse_regex_float_match,
    get_cell_dataframe_TS



include("StructsAndTypes.jl")
include("GetCellDataFunctions.jl")
include("GetDataPerTimeFunctions.jl")
include("GetTimeRelevantFunctions.jl")
include("GetCompetitionScoreFunctions.jl")
include("HelperFunctions.jl")
end











