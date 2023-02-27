function get_Nₜ_vec(data::DataFrame)
    return data.N
end


function test_compute_pop_ode_sol(t,N₀,r,β) 
    βr = (β ./ r)
    rt2 = (r .* t) ./ 2 
    Nₛₒₗ = (βr .+ (√N₀ .- βr) .* exp.(rt2)) .^2
    time_to_death = compute_time_extinction(N₀,β,r)
    
    if iszero(time_to_death[1].im)
        Nₛₒₗ[findall(x -> x > time_to_death[1].re,t)] .= 0
    elseif !iszero(time_to_death[1].im)
        Nₛₒₗ = Nₛₒₗ
    else
        # nothing
    end
    return Nₛₒₗ
end



"""
    Compute the N₀ given a Float64 values, in this case
    the Float64 is the time to death
"""
function N₀(tₐ)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    d = log(2)/tₐ
    gᵣ = log(2)/1.15
    return (2d/l*gᵣ)^2 *(pi*a)
end


"""
    Compute the N₀ given the death rate
"""
function N₀(d::DeathRate)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    gᵣ = log(2)/1.15
    return (2*d.d/(l*gᵣ))^2 *(π*a)
end

"""
    Compute the N₀ given the time to death
"""
function N₀(d::DeathRateTime)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    d = log(2)/d.d
    gᵣ = log(2)/1.15
    return (2d/l*gᵣ)^2 *(π*a)
end

"""
    Compute the N₀ given the apoptosis (or any other time to death)
    time
"""
function N₀(tₐ::ApoptosisTime)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    d = log(2)/tₐ.tₐ
    gᵣ = log(2)/1.15
    return (2d/l*gᵣ)^2 *(pi*a)
end

"""
    Compute the N₀ given the R₀
"""
function N₀(R₀::Radius)
    r = 0.5
    l = 3.0
    a = cₐ(r,l)
    
    return π*R₀.r^2/a
end

"""
    Compute the radius given the population size
"""
function R₀(N₀)
    r = 0.5
    l = 3.0
    a = cₐ(r,l)
    return sqrt(N₀*a/π)
end



"""
Get time series of each neighbour per type of cell
"""
function get_time_series_labels_diff_to_cell(::MachineState,cell_data)
    
    time_variants = get_time_vec_variants(MachineState(),cell_data)
    unique_time = time_variants.unique_time_vec

    pop_count_per_dt = [get_population_counts_per_time_step(cell_data,ut) for ut in unique_time]
    
    ts_neighbour_to_target = [i.count_neighbours_diff_to_target for i in pop_count_per_dt]
    ts_neighbour_to_attacker = [i.count_neighbours_diff_to_attacker for i in pop_count_per_dt]

    return CellLabelsTimeSeries(ts_neighbour_to_target,ts_neighbour_to_attacker)
end





"""
    Get cell population time series, as well as 
    get time series of neighbours different to cell type.
"""
function get_cell_population_ts(data,growth_rateₜ,growth_rateₐ)
    @chain data begin
        readdlm(_)
        get_cell_data(MachineState(),_)
        @aside @chain _ begin
            cell_ts = get_population_competition_time_series(
                MachineState(),_,
                growth_rateₜ,
                growth_rateₐ) 
        end
        @aside @chain _ begin
            cell_label_ts = get_time_series_labels_diff_to_cell(MachineState(),_)
        end
    end
    return CellTimeSeries(cell_ts,cell_label_ts)
end






"""
    Compute the count of neighbours to targets at a vector of time slice indices.
"""
function get_count_neighbour_targets(data::CellTimeSeries,time_slices)
    p⃗₀ = data.MachineStateTimeSeries.target    
    np⃗₀ = data.CellLabelsTimeSeries.neighbour_to_target ./ p⃗₀
    return np⃗₀[time_slices]
end


"""
    Take a data pather, filter to appropriate data paths,
    load data and get time series.
    Use `get_count_neighbour_targets` functions to get count of neighbours 
    to targets at time indices.
    Get mean values
"""
function get_mean_neighbours_at_time_time_slices(
    data_path,filter_pattern,
    growth_rateₜ,growth_rateₐ,time_slices)
    mean_neighbours = @chain data_path begin
        filter(x -> occursin(filter_pattern,x),_)    
        [get_cell_population_ts(i,growth_rateₜ,growth_rateₐ) 
        for i in _]
        [get_count_neighbour_targets(i,time_slices) for i in _]
        reduce(hcat,_)
        mean(_,dims=2)    
    end
    return mean_neighbours
end
    



"""
    Take data path with growth rates of both targets and attackers
    and returns the time series and older model of the competition.
"""
function simulation_population_TS(::MachineState,path,growth_rate₀,growth_rateₜ)
    Nₜ = @chain path begin
        readdlm(_) 
        get_cell_data(MachineState(),_) 
        get_population_competition_time_series(
                MachineState(),_,
                growth_rate₀,growth_rateₜ)
    end
    return Nₜ
end

"""
    Compute the time series of the targets population per 
    simulation and analytical solution.
"""
function simulation_analytical(
    path,
    growth_rate₀,
    growth_rateₜ,    
    radius,
    time_for_death)

    Ñ₀ = N₀(Radius(radius))

    β̃  = β(dᵣ(TimeForDeathRate(time_for_death)))
    Ñₜ = simulation_population_TS(
        MachineState(),path,growth_rate₀,growth_rateₜ)
    Nₛₒₗ = compute_pop_ode_sol(Ñₜ.time,Ñₜ.target[1],r,β̃ ) 
    
    
    time_to_death = compute_time_extinction.(Ñ₀,β̃ ,growth_rate₀)
    if iszero(time_to_death[1].im)
        Nₛₒₗ[findall(x -> x > time_to_death[1].re,Ñₜ.time)] .= 0
    elseif !iszero(time_to_death[1].im)
        Nₛₒₗ = Nₛₒₗ
    else
        # nothing
    end
    return (N_sim = Ñₜ,N_ana = Nₛₒₗ)
end








"""
    NOT GENERALISED
    Specific to path such that folder names separated by `/` or `\\` are removed by four times,then the remaining path removed all folder names until a single folder name remains. Here the radius to detect a neigbour is displayed.
"""
function get_neighbour_radius(data_path)
    neighbour_radius_str = @chain data_path begin
        dirname(_)
        dirname(_)
        dirname(_)
        dirname(_)
        basename(_)
    end
    return neighbour_radius_str
end


"""
    NOT GENERALISED
    Specific to path such that folder names separated by `/` or `\\` are removed by three times,then the remaining path removed all folder names until a single folder name remains. Here the radius to determine the size of the targets is displayed.
"""
function get_radius(data_path) 
    radius_str = @chain data_path begin
        dirname(_)
        dirname(_)
        dirname(_)
        basename(_)
        replace(_,"Radius_" => "")
        parse(Float64,_)
    end
    return radius_str
end





"""
    Get total population time series, simulated N
"""
function get_Nt(data::DataFrame)
    return @chain data begin
                @by([:cell_type_label,:time],
                    :N = size(:cell_type_label,1))
    end
end

"""
    Get target population from the the total population time series    
"""
function get_Nₜ(data::DataFrame)
    return @chain data begin
        @subset @byrow :cell_type_label == 0
    end
end



"""
    Get attacker population from the the total population time series
"""
function get_Nₐ(data::DataFrame)
    return @chain data begin
        @subset @byrow :cell_type_label == 1
    end
end

"""
    Compute the ODE solutions using the simulation data to grab the time duration.
"""
function compute_pop_ode_sol(
    Ñₜ::DataFrame,Δt,Ñ₀,β̃,r̃ₜ)
    t = range(minimum(Ñₜ.time),maximum(Ñₜ.time),step = Δt) 
    time_to_death = compute_time_extinction(Ñ₀,β̃,r̃ₜ)
    Nₛₒₗ = ChasteCellMachineReader.compute_pop_ode_sol(t,Ñ₀,r̃ₜ,β̃)
    
    if iszero(time_to_death[1].im)
        Nₛₒₗ[findall(x -> x > time_to_death[1].re,t)] .= 0
    elseif !iszero(time_to_death[1].im)
        Nₛₒₗ = Nₛₒₗ
    else
        # nothing
    end

    return DataFrame(
                time = t, 
                N = Nₛₒₗ)
end


"""
    Compute the N(t) ode of the per capita targets as they sit on 
    the boundary of the circle between the targets and attackers.
"""
function compute_pop_ode_sol(t,N₀,r,β) 
    βr = (β ./ r)
    rt2 = (r.*t) ./ 2 
    Nₛₒₗ = (βr .+ (√N₀ .- βr) .* exp.(rt2)) .^2
    time_to_death = compute_time_extinction(N₀,β,r)
    
    if iszero(time_to_death[1].im)
        Nₛₒₗ[findall(x -> x > time_to_death[1].re,t)] .= 0
    elseif !iszero(time_to_death[1].im)
        Nₛₒₗ = Nₛₒₗ
    else
        # nothing
    end
    return Nₛₒₗ
end



"""
    Compute the mean population with standard error bands.
    If one wants to obtain for targets or attackers, simply
    take that subset prior to calling this functions as it does
    not do so.
"""
function get_mean_pop_with_std_errors(data::DataFrame)
    return @chain data begin
                reduce(vcat,_.Nₜ) 
                @by(:time,
                :Ñ = mean(:N),
                :sdÑ = std(:N),
                :n = :realisations)
                @rtransform :sdÑ = isnan(:sdÑ) ? 0 : :sdÑ 
                @rtransform :Ñ₊ₛₑ = :Ñ + std_error(:sdÑ,:n)
                @rtransform :Ñ₋ₛₑ = :Ñ - std_error(:sdÑ,:n)
        end
end


"""
Compute the mean population with standard error bands.
If one wants to obtain for targets or attackers, simply
take that subset prior to calling this functions as it does
not do so. Using a dynamic input variable
"""
function get_mean_pop_with_std_errors(
    data)
    return @chain data begin
                reduce(vcat,_.Nₜ) 
                @by(:time,
                :Ñ = mean(:N),
                :sdÑ = std(:N),
                :n = :realisations)
                @rtransform :sdÑ = isnan(:sdÑ) ? 0 : :sdÑ 
                @rtransform :Ñ₊ₛₑ = :Ñ + std_error(:sdÑ,:n)
                @rtransform :Ñ₋ₛₑ = :Ñ - std_error(:sdÑ,:n)
        end
end


"""
    Set the initial population with the initial population
    Allows to keep the same interface
"""
function N₀(n₀::N₀AlreadySet)
    return n₀.n₀
end
