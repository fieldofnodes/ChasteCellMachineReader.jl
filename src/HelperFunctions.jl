

"""
    The area a capule with a cylinder length l and radius r
"""
cₐ(r,l) = π*r^2 + 2*r*l


"""
    The length a capule with a cylinder length l and radius r
"""
cₗ(r,l)  = 2*r + l


"""
    Using Dates package, get today's date and convert to YYYMMDD format
"""
today_nog_gap() = replace(today() |> string,"-" => "")



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
    Normalise a vector, x such that x∈[0,1]
"""
function get_norm_vec(x) 
    return (x .- minimum(x)) ./ (maximum(x) .- minimum(x))
end

"""
    Shift a U(0,1) to U(a,b)
"""
function shift_scale_uniform(min_value,max_value)
    range_dif = abs(min_value - max_value)
    return range_dif*rand()+min_value
end

"""
    Get mean values from a U(a,b)
"""
function get_mean_cycle_time(cycle_min_time,cycle_max_time,runs)
    return mean([shift_scale_uniform(cycle_min_time,cycle_max_time) for i in 1:runs])
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
    Take data::CellTimeSeries and plot competition per neighbours
"""
function plot_population_neighbours_competition(data::CellTimeSeries,radius)

    t⃗ = data.MachineStateTimeSeries.time       
    p⃗₀ = data.MachineStateTimeSeries.label₀
    p⃗₁ = data.MachineStateTimeSeries.label₁
    c⃗ = data.MachineStateTimeSeries.competition
    np⃗₀ = data.CellLabelsTimeSeries.neighbour_to_target ./ p⃗₀
    np⃗₁ = data.CellLabelsTimeSeries.neighbour_to_attacker ./ p⃗₁


    f = Figure(resolution = (800,600))

    ax1 = Axis(f[1,1],xlabel = "time",ylabel = "Absolute count", title = "1: Absolute Population Count")
    ax2 = Axis(f[1,2],xlabel = "time",ylabel = "Normalised count", title = "2: Count of Neighbours per Cell")
    ax3 = Axis(f[2,1:2],xlabel = "time", ylabel = "Death rate",title = "3: Death rate of targets")

    lp₀ = lines!(ax1,t⃗,p⃗₀)
    lp₁ = lines!(ax1,t⃗,p⃗₁)
    lnp₀ = lines!(ax2,t⃗,np⃗₁) # Count of targets per attackers
    lnp₁ = lines!(ax2,t⃗,np⃗₀) # Count of attackers per targets
    lc = lines!(ax3,t⃗,c⃗)
    
    Legend(
        f[3,1:2],
        [[lp₀,lnp₀],[lp₁,lnp₁],[lc]],
        ["Attackers","Targets","Death rate"],
        "Radius value: $(radius)",orientation = :horizontal)
    return f
end


"""
    Convert a matrix to a dataframe.
    Use the `:auto` variable.
    `:x1` is the radius of the target boundary
    `:x2`` is the competition value
    Then the variables from `:x3` onwards is shifted
    to `:x1` and so on.
"""
function convert_mat_to_df(mat)
    sz_mat = size(mat,2)
    df = @chain mat begin
        DataFrame(_,:auto) 
        rename(_,:x1 => :radius)
        rename(_,:x2 => :comp)
        rename(_,["x$i" => "x$(i-2)" for i in 3:sz_mat])
    end
    return df
end


"""
    Take a dataframe generated from `convert_mat_to_df()`.
    Generate a line of best fit using the lm model from the GLM 
    package.
"""
function compute_lm_predict(df)
    d̂f = []
    R² = []
    for x ∈ names(df[:,Not([:radius,:comp])])
        model_formula = term(:comp) ~ term((x))
        model = lm(model_formula,df)         
        ŷ = predict(model, df)
        push!(R²,r2(model))
        push!(d̂f,ŷ)
    end
    predict_df = DataFrame(d̂f,:auto)

    return (predict_df = predict_df,R² = R²)
end

"""
    Take two vectors, `x` and `y` and compute the line of best fit
    manually. Return the line of `ŷ`.
"""
function fit_line(x,y)
    ybar = mean(y)
    xbar = mean(x)
    xbar_diff = x .- xbar
    ybar_diff = y .- ybar
    slope = sum(xbar_diff .* ybar_diff)/sum(xbar_diff .^ 2)
    intercept = ybar - slope*xbar
    y_line = intercept .+ slope .* x
    return y_line
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
    NOT GENERALISED
    Specific to path such that folder names separated by `/` or `\\` are removed by two times,then the remaining path removed all folder names until a single folder name remains. Here the realisation number is displayed.
"""
function get_realisation_number(data_path)
    realisation = @chain data_path begin
        dirname(_)
        dirname(_)
        basename(_)
        replace(_,"Realisation_" => "")
        parse(Float64,_)
    end
    return realisation
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
    Take a dataframe and the number of realisations and includes this as a 
    column in the dataframe.
"""
function add_realisations(
    data::DataFrame,
    realisations::Int64)
    data.realisations .= realisations
    return data
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
    Compute the standard error from the standard deviation
    and the count of observartions.
"""
function std_error(σ,n)
    return σ/√n
end



"""
    If a value is missing, then output a zero
    Missing -> 0
"""
missing_to_zero(m::Missing) = 0


"""
    Function that takes an input and if missing than outputs 0 otherwise
    return what that number is.
"""
function if_missing_then_zero(x) 
    @assert isreal(x)|ismissing(x) "Input should be a Real number of missing"
    x̂ = ismissing(x) ? missing_to_zero(x) : x
    return x̂
end



"""
    Functions specifice to the chaste cell reader where we convert
    values per cell type label (0,1,3) and outputs to zero
    Function needs to be generalised.
"""
function convert_missing_to_zero(data) 
    df = @chain data begin 
            unstack(_,:cell_type_label,:N,renamecols=x->Symbol(:label_, x))
            @rtransform :label_0 = if_missing_then_zero(:label_0)
            @rtransform :label_3 = if_missing_then_zero(:label_3) 
            @rtransform :label_target = :label_0 + :label_3
            @select :time :label_1 :label_target
            stack(_,[:label_1, :label_target],:time)
            @rtransform :variable = :variable == "label_1" ? 1 : 0
            rename(_,:variable => :cell_type_label) 
            rename(_,:value => :N)
            @select :cell_type_label :time :N
        end
        return df
end




"""
    A once-off look behind patter when we want to change the patter
    of the look behind but still search for [0-9.]{1.5} after
"""
function input_neg_lookbehind(pattern)
    rx = Regex(string("(?<=",pattern,")[0-9.]{1,5}"))
    return rx
end

"""
    Take the regex patter, the string and the regex function to parse a 
    character into a Float64 type.
"""
function parse_regex_float_match(
    string_for_match,
    pattern_to_match,
    regex_match_function)

    rx = regex_match_function(pattern_to_match)
    m = match(rx,string_for_match,).match
    
    return parse(Float64,m)
end

"""
    Set the initial population with the initial population
    Allows to keep the same interface
"""
function N₀(n₀::N₀AlreadySet)
    return n₀.n₀
end
