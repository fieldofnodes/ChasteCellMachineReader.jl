####################################################################
# Title  : Load a range of simulations for a parameter exploration
# Author : Jonathan Miller         
# Date   : 20230304
# Aim    : Plot optimisations of data to parameter
#        : This script originally is fro CircularDomain
####################################################################

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
using Statistics

// # Functions
    """
        Take a dataframe which is a grouping and compute the mean for each column
        to return a dataframe of the same dim and a single dataframe.
    """
    function Statistics.mean(dat::Union{SubDataFrame,Any})
        data = dat.data
        names_dat = names(data[1])
        
        mat = []
        for n in eachindex(names_dat)
            push!(mat,mean([d[!,names_dat[n]] for d in data]))
        end
        return DataFrame(mat,names_dat)
    end


    """        
    plot_optimal_solution(
        sim_data,
        opt_data_local,opt_data_global,
        param_label,param_value)
    :sim_data: DataFrame - simulation time series
    :opt_data_local: NamedTuple - solution,β,d for global solution
    :opt_data_global: NamedTuple - solution,β,d for global solution
    :param_label: String - Parameter label
    :param_value: Float64 - parameter value

    Returns CairoMakie based plot

    """
    function plot_optimal_solution(
        sim_data,
        opt_data_local,opt_data_global,
        param_label,param_value)
        t = sim_data.time
        N = sim_data.target
        opt_sol_local = opt_data_local[:opt_sol]
        opt_param_local = opt_data_local[:β]
        d_β_local = opt_data_local[:d]
        opt_sol_global = opt_data_global[:opt_sol]
        opt_param_global = opt_data_global[:β]
        d_β_global = opt_data_global[:d]
    
        f = Figure(resolution = (800,600))
            ax = Axis(f[1,1],xlabel = L"T", ylabel = L"N")
            lines!(ax,t,N,label = "data, $(param_label) = $(param_value)")
            lines!(ax,t,opt_sol_local, label="Local βₒₚₜ ≈ $(round(opt_param_local,digits=3))\nd_from_β = $(round(d_β_local,digits=3))")
            lines!(ax,t,opt_sol_global, label="Global βₒₚₜ ≈ $(round(opt_param_global,digits=3))\nd_from_β = $(round(d_β_global,digits=3))")
            axislegend()
        return f
    end

    function save_sim_TS_plots(
        plots,
        plot_path,
        sim_type,
        param_label,
        param_value)

        each_plot = "plot_".*[
            "cell_populations",
            "machines_diff_state",
            "per_cap_machines_state",
            "total_machines",
            "per_cap_total_machines",
            "total_neighbours",
            "per_cap_total_neighbours"]
    
        file_ext = "png"

        folder_dir = string(plot_path,param_label,"_",param_value,"/")
        
        if !isdir(folder_dir)
            mkdir(folder_dir)
        end
        
        output_path = string.(
            folder_dir,
            today_nog_gap(),"_",
            each_plot,"_",
            sim_type,"_",
            param_label,"_",
            param_value,".",
            file_ext)

        save.(output_path,plots)
        return "Saved to file"
    end

    function save_sim_fitted_plots(
        opt_local,
        opt_global,
        fitted_plots,
        plot_path,
        sim_type,
        param_label,
        param_value)

        folder_dir = string(plot_path,param_label,"_",param_value,"/")
        
        if !isdir(folder_dir)
            mkdir(folder_dir)
        end
   
        opt_param_local = opt_local[:β]
        opt_param_global = opt_global[:β]

        filename = 
            string(
                folder_dir,
                today_nog_gap(),"_",
                sim_type,"_",
                "Local-Global_",
                param_label,
                "-",
                param_value,
                "_beta-local-",
                round(opt_param_local,digits=3),
                "-global-",
                round(opt_param_global,digits=3),
                "_mean_TS.png")

        save(filename,fitted_plots)

        return "Saved to file"
    end

    function save_fitted_range_plots(
        plots,
        plot_path,
        param_label,
        sim_type
        )

        folder_dir = string(plot_path,param_label,"_","Range","/")
            
        if !isdir(folder_dir)
            mkdir(folder_dir)
        end

        
        
        each_plot = [
            "population_fitted_opt",
            "population_fitted_beta",
            "population_fitted_death_rate"]
            
        file_ext = "png"

        output_path = string.(
            folder_dir,
            today_nog_gap(),"_",
            each_plot,"_",
            sim_type,".",
            file_ext)


        save.(output_path,plots)
        
        return "Saved to file"
    end

    function write_data_to_folder(
        data,
        cell_path,
        sim_type,
        param_label,
        param_value)
        
        new_dir = @chain cell_path begin
            dirname(_)
            dirname(_) 
            replace(_,"chaste_" => "")
        end
    
        if !isdir(dirname(new_dir))
            mkdir(dirname(new_dir))
        end

        if !isdir(new_dir)
            mkdir(new_dir)
        end
        
        output_desc = "TS_df"
        file_ext = "csv"
    
        output_path = string(
            new_dir,"/",
            today_nog_gap(),"_",
            output_desc,"_",
            sim_type,"_",
            param_label,"_",
            param_value,".",
            file_ext)
    
        CSV.write(output_path,data)
        return "Saved to file"
    
    end
    
    
//

// # Input parameters
    sim_type = "CircularDomain"
    chaste_data = "chaste_data/"
    data_path = "data/"
    plot_path = "plots/"
    find_pattern = "machinestate.dat"
    param_label = "K1"
    param_label_in_chase = "K_1"
    rx_find_pattern = r"K_1_[0-9.]{1,5}"
    rx_param_label_pattern = r"K_1"
    rx_param_value_pattern = r"([0-9]{2,3}\.[0-9]{1})|([0-9]{1,3}$)"
    rx_realisation_pattern = r"Realisation_[0-9]{1,2}"
    rx_ralisation_label = r"[0-9]{1,2}$"
    Δt = 1.0/1200.0
    target_time_for_death = 0.25
    ϵ = 0.2
    t_t6ss = 1/20
    cycle_time_01 = 0.5
    cycle_time_02 = 1.5
    target_mean_cycle_time = (cycle_time_01 + cycle_time_02)/2
    p = [target_time_for_death,target_mean_cycle_time]
    
    
//



// # Cell paths, data, optimisation and save data to file
    """
        Get cell paths, data and optimisation
    """
    cell_path = @chain chaste_data begin
        find(_,0,7,find_pattern) 
    end

    """
    Process folder names for saving data appropriateley
    Manually enter the specific regex we want to extract values for
    """
    @chain cell_path begin
        @aside @chain _ begin
            [match(rx_find_pattern,m).match for m in _] 
            @aside @chain _ begin
                label = [replace(match.(rx_param_label_pattern,i).match,"_"=>"") for i in _]
            end
            @aside @chain _ begin
                value = [parse(Float64,
                        match(rx_param_value_pattern,i).match) 
                        for i in _]
            end
        end 
        [match(rx_realisation_pattern,m).match
                for m in _]
        @aside @chain _ begin
            realisation = [parse(Float64,
                match(rx_ralisation_label,m).match) for m in _]
        end
    end

    
    """
        Compute metrics based on data and parametes.
    """
    death_rate = 
        [get_death_rate(K₁(p),
            target_time_for_death,
            t_t6ss,
            target_mean_cycle_time,ϵ).d for p in value]

    β⃗ = β.(DeathRate.(death_rate))
    cdfs = get_cell_dataframe_TS.(cell_path)
    opt_local = [get_opt_param_sol(Local(),c,p) for c in cdfs]
    opt_global = [get_opt_param_sol(Global(),c,p) for c in cdfs]
    
    
    """
        Plot mean population results 
        For each time series caculate the associate mean time series per
        parameters 
        Plot the normal plots for each param
        Plot the optimisation sol for each parameter
        Combine into a single dataframe for each parameter

    """    
    mn_df = @chain cdfs begin
        DataFrame(param_value = value,data = _)
        groupby(_,:param_value)
        [mean(p) for p in _]
        DataFrame(param_value = unique(value), mean_dfs = _)
        @orderby :param_value
        @rtransform :opt_local = get_opt_param_sol(Local(),:mean_dfs,p)
        @rtransform :opt_global = get_opt_param_sol(Global(),:mean_dfs,p)
        @rtransform :param_label = param_label
        @rtransform :plot = plot_cell_populations_TS(:mean_dfs)
        @rtransform :fitted_plot = plot_optimal_solution(
                            :mean_dfs,
                            :opt_local,
                            :opt_global,
                            :param_label,
                            :param_value)
        @rtransform :save_plot = save_sim_TS_plots(
                            :plot,
                            plot_path,
                            sim_type,
                            :param_label,
                            :param_value)
        @rtransform :save_fitted_plot = save_sim_fitted_plots(
                            :opt_local,
                            :opt_global,
                            :fitted_plot,
                            plot_path,
                            sim_type,
                            :param_label,
                            :param_value)
    end

    
    


    """
        Aggregate data
        Gather all solutions into a single data frame.
        We collect the single parameter, population values, β and death rates.
    """
    df = @chain begin
        DataFrame(
            param_value = value,
            NtN0 = [i.target[end]/i.target[1] for i in cdfs],
            opt_NtN0_local = [o.opt_sol[end]/o.opt_sol[1] for o in opt_local],
            opt_NtN0_global = [o.opt_sol[end]/o.opt_sol[1] for o in opt_global],
            β_local = [i.β for i in opt_local],
            β_global = [i.β for i in opt_global],
            β⃗_from_d = β⃗,
            d_local = [i.d for i in opt_local],
            d_global = [i.d for i in opt_global],
            theory_d = death_rate,)
        @orderby :param_value
    end
    
    """
        We order by the parameter, then group by it and 
        compute the mean for each remaining variable.
    """
    gdf = @chain df begin
        @by(:param_value,
            :NtN0 = mean(:NtN0),
            :opt_NtN0_local = mean(:opt_NtN0_local),
            :opt_NtN0_global = mean(:opt_NtN0_global),
            :β_local = mean(:β_local),
            :β_global = mean(:β_global),
            :β⃗_from_d = mean(:β⃗_from_d),
            :d_local = mean(:d_local),
            :d_global = mean(:d_global),
            :theory_d = mean(:theory_d))
    end
      

    """
    Generate and write to file each TS per chaste file
    We will save each data frame, plus the two aggregated dataframes
    :cdfs,mn_df,df,gdf -TODO
    """


    """
        Save all cell data to file
    """
    [write_data_to_folder(
        cdfs[i],
        cell_path[i],
        sim_type,
        param_label,
        value[i]) for i in eachindex(cdfs)]



    @chain mn_df begin
        @select :param_label :param_value :mean_dfs
        @rtransform :param_label = replace(:param_label, "K1" => "K_1")
        @rtransform :param_label_value = string(
                :param_label,"_",
                :param_value)
        @rtransform :param_label_value = replace(
                :param_label_value,
                "K_1_1.0" =>"K_1_1")
        @rtransform :param_label_value = replace(
                :param_label_value,
                "K_1_100.0" =>"K_1_100") 
        @rtransform :path = string(
                data_path,
                :param_label_value,"/",
                today_nog_gap(),"_",
                "TS_mean_values_df","_",
                :param_label_value,"_",
                sim_type,
                ".csv")
        @select :path :mean_dfs
        @rtransform :saved = CSV.write(:path,:mean_dfs)
    end
        
    """
    For df all observations with the param values and fitted
    Mean is taken from these    
    """
    path = string(data_path,today_nog_gap(),"_",
                "param_fitted_all_K1_obvs.csv")
    
    CSV.write(path,df)

    """
    For df all observations with the param values and fitted
    Mean is taken from these    
    """
    path = string(data_path,today_nog_gap(),"_",
                "param_fitted_per_each_K1_group_obvs.csv")
    
    CSV.write(path,gdf)


        

    

//




// # Generate plotted results and save to file

    """
        Plot Population at time, t, divided by population at time, 0, (Nt/N0) as a function of K1
        Do this for the data, the local and the global optimisations
    """
    pop_fig = Figure(resolution = (800,600))
        ax = Axis(pop_fig[1,1],xlabel = param_label, ylabel = L"N_t/N_0")
        param = gdf.param_value
        NtN0 = gdf.NtN0
        oNtN0_l = gdf.opt_NtN0_local
        oNtN0_g = gdf.opt_NtN0_global
        scatter!(ax,param,NtN0,label="Simulation")
        scatter!(ax,param,oNtN0_l,label="Local Opt")
        scatter!(ax,param,oNtN0_g,label="Global Opt")
        axislegend()
    pop_fig


    """
        Plot β as a function of K1 
        Do this for the local and the global
    """
    β_fig = Figure(resolution = (800,600))
        ax = Axis(β_fig[1,1],xlabel = param_label, ylabel = L"\beta")
        param = gdf.param_value
        β_local = gdf.β_local
        β_global = gdf.β_global
        β⃗_from_d = gdf.β⃗_from_d
        scatter!(ax,param,β⃗_from_d,label="β Theory")
        scatter!(ax,param,β_local,label="Local Opt")
        scatter!(ax,param,β_global,label="Global Opt")
        lines!(ax,param,[maximum(β⃗)])
        axislegend(position=:rb)
    β_fig


    """
        Plot d as a function of K1 
        Do this for the local and the global
    """
    d_fig = Figure(resolution = (800,600))
        ax = Axis(d_fig[1,1],xlabel = param_label, ylabel = L"d")
        param = gdf.param_value
        d_local = gdf.d_local
        d_global = gdf.d_global
        theory_d = gdf.theory_d
        scatter!(ax,param,theory_d,label="D Theory")
        scatter!(ax,param,d_local,label="Local Opt")
        scatter!(ax,param,d_global,label="Global Opt")
        lines!(ax,param,[log(2)/0.25])
        axislegend(position=:rb)
    d_fig

    plots = [pop_fig, β_fig, d_fig]

    save_fitted_range_plots(
        plots,
        plot_path,
        param_label,
        sim_type
        )





//

 
