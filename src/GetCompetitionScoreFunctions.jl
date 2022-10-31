"""
    Take two time series and a time vector to compute the competition score. This is score is normalised according to the initial population count per each population.

    The quation for this score is defined as a time series,

    ```math c(t) = \frac{1}{t}\log\left( \frac{s(t)_{2}}{s(t=0)_{2}}\frac{s(t)_{1}}{s(t=0)_{1}} \right) ``

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
    Once the cell data has been read into the vector of `MachineStateCellProperties`, then this function computes the time vector which adjusts the initial time value to `0`. 
    
    Each cell label is extracted and the according time series of total population counts is defined. 

    The function returns the time shifted time series, the time series for cell label `0` and `1` and the time series competition.
"""
function get_population_competition_time_series(::MachineState,cell_data)
    time_vars = get_time_vec_variants(MachineState(),cell_data)

    population_counts_per_time_step = [get_population_counts_per_time_step(cell_data,t) for t in time_vars.unique_time_vec]
    time_series₁ = [p.count_label_0 for p in population_counts_per_time_step]
    time_series₂ = [p.count_label_1 for p in population_counts_per_time_step]
    competition = get_competition_score_function(time_series₁,time_series₂,time_vars.unique_shifted_time_vec)



    return MachineStateTimeSeries(time_vars.unique_shifted_time_vec,time_series₁,time_series₂,competition)
end
