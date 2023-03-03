"""
    Sample script to take a folder of chaste data
    create a data frame of the time series for each simulation
    create some plots
    save to file
"""

using Pkg
Pkg.activate(".")
using ChasteCellMachineReader
using Find
using DataFrames
using DataFramesMeta
using CSV
using CairoMakie

// # Input parameters
    sim_type = "CircularDomain"
    chaste_data = "chaste_data/"
    data_path = "data/"
    plot_path = "plots/"
    find_pattern = "machinestate.dat"
    
    Δt = 1.0/1200.0
    target_time_for_death = 0.25
    ϵ = 0.2
    t_t6ss = 1/20
    cycle_time_01 = 0.5
    cycle_time_02 = 1.5
    target_mean_cycle_time = (cycle_time_01 + cycle_time_02)/2
    p = [target_time_for_death,target_mean_cycle_time]
    param = (label = "K1", value = 100)
    

    
    death_rate = 
        get_death_rate(K₁(param[:value]),
            target_time_for_death,
            t_t6ss,
            target_mean_cycle_time,ϵ)

    d = round(death_rate.d,digits=3)
    β̂  = round(β(death_rate),digits = 3)
    
    
//



// # Cell paths, data, optimisation and save data to file
    """
        Get cell paths, data and optimisation
    """
    cell_path = @chain chaste_data begin
        find(_,0,7,find_pattern) 
    end

    cdfs = get_cell_dataframe_TS.(cell_path[1])
    opt_local = get_opt_param_sol(Local(),cdfs,p)
    opt_global = get_opt_param_sol(Global(),cdfs,p)

    """
    Generate and write to file each TS per chaste file
    """
    output_desc = "TS_df"
    file_ext = "csv"

    output_path = string.(
        data_path,
        today_nog_gap(),"_",
        output_desc,"_",
        sim_type,"_",
        param[:label],"_",
        param[:value],".",
        file_ext)

    CSV.write(output_path,cdfs)


//




// # Generate plotted results and save to file

    """
        Save plots to file
    """
    plots = plot_cell_populations_TS(cdfs)

    each_plot = "plot_".*[
        "cell_populations",
        "machines_diff_state",
        "per_cap_machines_state",
        "total_machines",
        "per_cap_total_machines",
        "total_neighbours",
        "per_cap_total_neighbours"]

    file_ext = "png"

    output_path = string.(
        plot_path,
        today_nog_gap(),"_",
        each_plot,"_",
        sim_type,"_",
        param[:label],"_",
        param[:value],".",
        file_ext)


    save.(output_path,plots)

    
    
    f = Figure(resolution = (800,600))
        t = cdfs.time
        N = cdfs.target
        opt_sol_local = opt_local[:opt_sol]
        opt_param_local = opt_local[:β]
        d_β_local = opt_local[:d]
        opt_sol_global = opt_global[:opt_sol]
        opt_param_global = opt_global[:β]
        d_β_global = opt_global[:d]
        ax = Axis(f[1,1],xlabel = L"T", ylabel = L"N")
        lines!(ax,t,N,label = "data, $(param[:label]) = $(param[:value])")
        lines!(ax,t,opt_sol_local, label="Local βₒₚₜ ≈ $(round(opt_param_local,digits=3))\nd_from_β = $(round(d_β_local,digits=3))")
        lines!(ax,t,opt_sol_global, label="Global βₒₚₜ ≈ $(round(opt_param_global,digits=3))\nd_from_β = $(round(d_β_global,digits=3))")
        axislegend("Inference\nd from $(param[:label]) = $(d)\nβ from d = $(β̂ )")
    f

    filename = 
        string(
            "plots/",
            today_nog_gap(),"_",
            "Local-Global_",
            param[:label],
            "-",
            param[:value],
            "_beta-local-",
            round(opt_param_local,digits=3),
            "-global-",
            round(opt_param_global,digits=3),
            "_mean_TS.png")

    save(filename,f)
    

//

 
