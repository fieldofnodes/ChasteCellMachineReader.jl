"""
    Compute a time series and shift every element such that the initial one is `0`.
"""
function shift_time_back_to_zero(time)
    return time .- time[1]
end


"""
    Compute different expressions of the time series
    1. `time_vec`
    2. `unique_time_vec`
    3. `shifted_time_vec`
    4. `unique_shifted_time_vec`
"""
function get_time_vec_variants(::MachineState,cell_data)
    time_vec = [c.time for c in cell_data]
    shifted_time_vec = shift_time_back_to_zero(time_vec)
    unique_shifted_time_vec = shifted_time_vec |> unique
    unique_time_vec = time_vec |> unique
    
    return TimeVecVariants( 
        time_vec,
        unique_time_vec,
        shifted_time_vec,
        unique_shifted_time_vec)
end



function get_time_vec(data::DataFrame)
    return data.time
end
