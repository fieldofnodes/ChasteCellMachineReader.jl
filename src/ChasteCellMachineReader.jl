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
@reexport using Optimization
@reexport using OptimizationOptimJL
@reexport using OptimizationBBO


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
    CellLabelsTimeSeries,
    CellTimeSeries,
    TimeForDeathRate,
    DeathRate,
    DeathRateTime,
    N₀AlreadySet,
    K₁,
    T6SSKineticsPositionTime,
    CellDivisionTime,
    Local,
    Global,
    DeathRateMichaelisMentenForm,
    CellBoundaryNormal,
    DataAllCellLabels,
    get_cell_data,
    get_cell_dataframe,
    get_cell_dataframe_TS,
    cₐ,
    cₗ,
    get_competition_score_function,
    get_population_competition_time_series,
    get_death_rate,
    get_time_to_die_β,
    dᵣ,
    β,
    get_mean_competition,
    compute_target_death_rate,
    get_competition_per_params,
    get_mean_comp,
    compute_time_extinction,
    get_type_count_per_time,
    get_population_counts_per_time_step,
    shift_time_back_to_zero,
    get_time_vec_variants,
    get_time_vec,
    today_nog_gap,
    get_norm_vec,
    shift_scale_uniform,
    get_mean_cycle_time,
    convert_mat_to_df,
    get_realisation_number,
    add_realisations,
    std_error,
    missing_to_zero,
    if_missing_then_zero,
    convert_missing_to_zero,
    input_neg_lookbehind,
    parse_regex_float_match,
    new_folder_name,
    generate_cell_output_path,
    make_output_path,
    convert_vecs_to_mat,
    optimise_β,
    optimise_K,
    first_analytical_solution,
    myloss,
    compute_lm_predict,
    fit_line,
    plot_cell_populations_TS,
    plot_population_neighbours_competition,
    get_Nₜ_vec,
    test_compute_pop_ode_sol,
    N₀,
    R₀,
    get_time_series_labels_diff_to_cell,
    get_cell_population_ts,
    get_count_neighbour_targets,
    get_mean_neighbours_at_time_time_slices,
    simulation_population_TS,
    simulation_analytical,
    get_neighbour_radius,
    get_radius,
    get_Nt,
    get_Nₜ,
    get_Nₐ,
    compute_pop_ode_sol,
    get_mean_pop_with_std_errors,
    get_opt_param_sol,
    extract_dataframe_names


include("StructsAndTypes.jl")
include("GetCellDataFunctions.jl")
include("GetCompetitionScoreFunctions.jl")
include("GetDataPerTimeFunctions.jl")
include("GetTimeRelevantFunctions.jl")
include("HelperFunctions.jl")
include("OptimisationFunctions.jl")
include("PlotCellData.jl")
include("PlottingFunctions.jl")
include("PopulationRelatedFunctions.jl")


end











