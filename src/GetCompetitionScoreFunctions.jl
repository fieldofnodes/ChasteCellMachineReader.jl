"""
    Take two time series and a time vector to compute the competition score. This is score is normalised according to the initial population count per each population.

    
    This function returns the finte competition time series.

"""
function get_competition_score_function(time_series₁,time_series₂,time)
    logs_vec = Float64.(
        log.(
            (time_series₂ ./ time_series₂[1]) ./ 
            (time_series₁ ./ time_series₁[1])) ./ 
            time)
    
    return logs_vec
end


"""
    Take two time series and a time vector to compute the competition score. This is score is normalised according to the initial population count per each population.

    
    This function returns the finte competition time series.

"""
function get_competition_score_function(::MachineState,time_series₁,time_series₂,time,growth_rate₁,growth_rate₂)
    basic_comp = Float64.(log.(
            (time_series₁[1] .* time_series₂) ./ (time_series₂[1] .* time_series₁)))

    comp = 1 .+ ((1 ./ time) .* ((1/growth_rate₁) .* basic_comp)) .- (growth_rate₁ / growth_rate₂)
    
    return comp
end




"""
    Once the cell data has been read into the vector of `MachineStateCellProperties`, then this function computes the time vector which adjusts the initial time value to `0`. 
    
    Each cell label is extracted and the according time series of total population counts is defined. 

    The function returns the time shifted time series, the time series for cell label `0` and `1` and the time series competition.
"""
function get_population_competition_time_series(::MachineState,cell_data,growth_rate₁,growth_rate₂)
    time_vars = get_time_vec_variants(MachineState(),cell_data)

    population_counts_per_time_step = [get_population_counts_per_time_step(cell_data,t) for t in time_vars.unique_time_vec]
    time_series₁ = [p.count_label_0 for p in population_counts_per_time_step]
    time_series₂ = [p.count_label_1 for p in population_counts_per_time_step]
   # competition = 
   #     get_competition_score_function(
   #         MachineState(),time_series₁,time_series₂,time_vars.unique_shifted_time_vec,growth_rate₁,growth_rate₂)



    return MachineStateTimeSeries(time_vars.unique_shifted_time_vec,time_series₁,time_series₂)
end


function get_death_rate(β)
    cᵣ = 0.5
    cylₗ = 3.0
    l = cₗ(cᵣ,cylₗ)
    a = cₐ(cᵣ,cylₗ)
    return DeathRate((β*l)/(2*√(π*a)))
end 

function get_death_rate(::CellBoundaryNormal,β)
    cᵣ = 0.5
    cylₗ = 3.0
    l = cᵣ*2
    a = cₐ(cᵣ,cylₗ)
    return DeathRate((β*l)/(2*√(π*a)))

end 


function get_time_to_die_β(d::DeathRate)
    return DeathRateTime(log(2)/d.d)
end     

function get_time_to_die_β(β)
    cᵣ = 0.5
    cylₗ = 3.0
    l = cₗ(cᵣ,cylₗ)
    a = cₐ(cᵣ,cylₗ)
    dₜ = (log(2)*2*√(π*a))/(β*l)
    return dₜ
end 



function get_death_rate(
    k₁::K₁,
    tₐ::ApoptosisTime,
    t_t6ss::T6SSKineticsPositionTime,
    t_cell_div::CellDivisionTime,
    ϵ)

    p₁ = π*√(2/ϵ)
    p₂ = t_t6ss.t/t_cell_div.t
    p₃ = log(2)/k₁.k₁
    return DeathRate(log(2)/(tₐ.tₐ + p₁*p₂*p₃))
end


function get_death_rate(
    ::DeathRateMichaelisMentenForm,
    k₁::Vector{K₁},
    tₐ::ApoptosisTime)
    d(K) = [log(2) ./ tₐ.tₐ .* (k.k₁ ./ (K .+ k.k₁)) for k in k₁]
    return d
end



function get_death_rate(
            k₁::K₁,
            target_time_for_death,
            t_t6ss,
            target_mean_cycle_time,ϵ)
    d = get_death_rate(k₁,
        ApoptosisTime(target_time_for_death),
        T6SSKineticsPositionTime(t_t6ss),
        CellDivisionTime(target_mean_cycle_time),
        ϵ)  
        return d
end

"""
    Compute the death rate given the N₀, ther growth rate and the 
    length and area of the cell.
"""
function dᵣ(N₀)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    gᵣ = log(2)/1.15
    return ((l*gᵣ)*√(N₀/(π*a)))/2
end


"""
    Compute the death rate from the time to death type
"""
function dᵣ(t::TimeForDeathRate)
    return DeathRate(log(2)/t.t)
end
    
"""
    Compute β from the death rate
"""
function β(d::DeathRate)
    cᵣ = 0.5
    cylₗ = 3.0
    l = cₗ(cᵣ,cylₗ)
    a = cₐ(cᵣ,cylₗ)
    return (2*d.d*√(π*a))/(l)
end



"""
    Takes data path to a folder of realisations.
    Filters data paths according to a filter_pattern.
    Computes time series for cell Population.
    Computed the so-called competition of the attackers and 
    takes the mean of the last element per time series.
    This assumes there is convergence.
"""
function get_mean_competition(
    data_path,filter_pattern,
    growth_rateₜ,growth_rateₐ)
    
    c̄ =@chain data_path begin
        filter(x -> occursin(filter_pattern,x),_)
        [get_cell_population_ts.(i,growth_rateₜ,growth_rateₐ).MachineStateTimeSeries.competition[end] for i in _]
        mean(_)
    end

    return c̄
end


"""
    Compute per capita death rate for targets on the circular
    boundary dividing the targets and attackers.
"""
function compute_target_death_rate(t,N,cell_cycle_time,l,a)
    r = log(2) ./ cell_cycle_time
    prod₁ = (l .* r) ./ (2 .* sqrt.(π .* a))
    prod₂ = (sqrt.(N) .- sqrt.(N[1]) .* exp.((r .* t) ./ 2)) ./ (1 .- exp.((r .* t) ./ 2))
    d = prod₁ .* prod₂
    return d
end


"""
    Load data path.
    Get cell data.
    Get the population time series and 
    competition time series.
"""
function get_competition_per_params(path,growth_rateₜ,growth_rateₐ)
    mat_data = readdlm(path)
    cell_data = get_cell_data(MachineState(),mat_data)
    time_series = 
        get_population_competition_time_series(
            MachineState(),
            cell_data,
            growth_rateₜ,growth_rateₐ)
    return time_series.competition[end]
end

"""
    Take a competition score vector.
    Compute the sum divided by the length of the vector.
    Return the arithmetic mean of the compeition score.
"""
function get_mean_comp(competition_score)
    return sum(competition_score)/length(competition_score)
end



"""
    Compute the time to extinction.
"""
function compute_time_extinction(N₀,β,r)
    twodivr = 2 ./ r
    βdivr = β ./ r
    divβr = -βdivr ./ (√N₀ .- βdivr )
    return twodivr .* log.(Complex.(divβr))
end
