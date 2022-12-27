

"""
    The area a capule with a cylinder length l and radius r
"""
cₐ(r,l) = π*r^2 + 2*r*l


"""
    The length a capule with a cylinder length l and radius r
"""
cₗ(r,l)  = 2*r + l



function N₀(tₐ)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    d = log(2)/tₐ
    gᵣ = log(2)/1.15
    return (2d/l*gᵣ)^2 *(pi*a)
end

function N₀(tₐ::ApoptosisTime)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    d = log(2)/tₐ.tₐ
    gᵣ = log(2)/1.15
    return (2d/l*gᵣ)^2 *(pi*a)
end

function N₀(R₀::Radius)
    r = 0.5
    l = 3.0
    a = cₐ(r,l)
    
    return π*R₀.r^2/a
end

function R₀(N₀)
    r = 0.5
    l = 3.0
    a = cₐ(r,l)
    return sqrt(N₀*a/π)
end



function get_norm_vec(x) 
    return (x .- minimum(x)) ./ (maximum(x) .- minimum(x))
end

function shift_scale_uniform(min_value,max_value)
    range_dif = abs(min_value - max_value)
    return range_dif*rand()+min_value
end


function get_mean_cycle_time(cycle_min_time,cycle_max_time,runs)
    return mean([shift_scale_uniform(cycle_min_time,cycle_max_time) for i in 1:runs])
end





function get_competition_per_params(path,growth_rate₁,growth_rate₂)
    mat_data = readdlm(path)
    cell_data = get_cell_data(MachineState(),mat_data)
    time_series = 
        get_population_competition_time_series(
            MachineState(),
            cell_data,
            growth_rate₁,growth_rate₂)
    return time_series.competition[end]
end


function get_mean_comp(competition_score)
    return sum(competition_score)/length(competition_score)
end

"""
Get time series of each neighbour per type of cell
"""
function get_time_series_labels_diff_to_cell(::MachineState,cell_data)
    
    time_variants = get_time_vec_variants(MachineState(),cell_data)
    unique_time = time_variants.unique_time_vec

    pop_count_per_dt = [get_population_counts_per_time_step(cell_data,ut) for ut in unique_time]
    
    ts_neighbour_to_label_0 = [i.count_neighbours_diff_to_cell_0 for i in pop_count_per_dt]
    ts_neighbour_to_label_1 = [i.count_neighbours_diff_to_cell_1 for i in pop_count_per_dt]

    return CellLabelsTimeSeries(ts_neighbour_to_label_0,ts_neighbour_to_label_1)
end


today_nog_gap() = replace(today() |> string,"-" => "")



function get_cell_population_ts(data,growth_rate₁,growth_rate₂)
    @chain data begin
        readdlm(_)
        get_cell_data(MachineState(),_)
        @aside @chain _ begin
            cell_ts = get_population_competition_time_series(
                MachineState(),_,
                growth_rate₁,
                growth_rate₂) 
        end
        @aside @chain _ begin
            cell_label_ts = get_time_series_labels_diff_to_cell(MachineState(),_)
        end
    end
    return CellTimeSeries(cell_ts,cell_label_ts)
end



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


function get_mean_competition(
    data_path,filter_pattern,
    growth_rate₁,growth_rate₂)
    
    c̄ =@chain data_path begin
        filter(x -> occursin(filter_pattern,x),_)
        [get_cell_population_ts.(i,growth_rate₁,growth_rate₂).MachineStateTimeSeries.competition[end] for i in _]
        mean(_)
    end

    return c̄
end

function get_count_neighbour_targets(data::CellTimeSeries,time_slices)
    p⃗₀ = data.MachineStateTimeSeries.label₀    
    np⃗₀ = data.CellLabelsTimeSeries.neighbour_to_label_0 ./ p⃗₀
    return np⃗₀[time_slices]
end


function get_mean_neighbours_at_time_time_slices(
    data_path,filter_pattern,
    growth_rate₁,growth_rate₂,time_slices)
    mean_neighbours = @chain data_path begin
        filter(x -> occursin(filter_pattern,x),_)    
        [get_cell_population_ts(i,growth_rate₁,growth_rate₂) 
        for i in _]
        [get_count_neighbour_targets(i,time_slices) for i in _]
        reduce(hcat,_)
        mean(_,dims=2)    
    end
    return mean_neighbours
end
    
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

function compute_target_death_rate(t,N,cell_cycle_time,l,a)
    r = log(2) ./ cell_cycle_time
    prod₁ = (l .* r) ./ (2 .* sqrt.(π .* a))
    prod₂ = (sqrt.(N) .- sqrt.(N[1]) .* exp.((r .* t) ./ 2)) ./ (1 .- exp.((r .* t) ./ 2))
    d = prod₁ .* prod₂
    return d
end




function plot_population_neighbours_competition(data::CellTimeSeries,radius)

    t⃗ = data.MachineStateTimeSeries.time       
    p⃗₀ = data.MachineStateTimeSeries.label₀
    p⃗₁ = data.MachineStateTimeSeries.label₁
    c⃗ = data.MachineStateTimeSeries.competition
    np⃗₀ = data.CellLabelsTimeSeries.neighbour_to_label_0 ./ p⃗₀
    np⃗₁ = data.CellLabelsTimeSeries.neighbour_to_label_1 ./ p⃗₁


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



function compute_pop_ode_sol(t,N₀,r,b) 
    br = (b ./ r)
    rt2 = (r.*t) ./ 2 
    pop = (br .+ (√N₀ .- br) .* exp.(rt2)) .^2
    return pop 
end

cₐ(r,l) = π*r^2 + 2*r*l
cₗ(r,l)  = 2*r + l

function N₀(d::DeathRate)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    gᵣ = log(2)/1.15
    return (2*d.d/(l*gᵣ))^2 *(π*a)
end


function N₀(d::DeathRateTime)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    d = log(2)/d.d
    gᵣ = log(2)/1.15
    return (2d/l*gᵣ)^2 *(π*a)
end


function N₀(tₐ::ApoptosisTime)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    d = log(2)/tₐ.tₐ
    gᵣ = log(2)/1.15
    return (2d/l*gᵣ)^2 *(pi*a)
end

function N₀(R₀::Radius)
    r = 0.5
    l = 3.0
    a = cₐ(r,l)
    
    return π*R₀.r^2/a
end

function R₀(N₀)
    r = 0.5
    l = 3.0
    a = cₐ(r,l)
    return sqrt(N₀*a/π)
end

function dᵣ(N₀)
    r = 0.5
    l = 3.0
    l = cₗ(r,l)
    a = cₐ(r,l)
    gᵣ = log(2)/1.15
    return ((l*gᵣ)*√(N₀/(π*a)))/2
end

function dᵣ(t::TimeForDeathRate)
    return DeathRate(log(2)/t.t)
end
    

function β(d::DeathRate)
    cᵣ = 0.5
    cylₗ = 3.0
    l = cₗ(cᵣ,cylₗ)
    a = cₐ(cᵣ,cylₗ)
    return (2*d.d*√(π*a))/(l)
end



function compute_time_extinction(N₀,β,r)
    twodivr = 2/r
    βdivr = β/r
    divβr = -βdivr/(√N₀ - βdivr )
    return twodivr * log(Complex(divβr))
end


function simulation_population_TS(::MachineState,path,growth_rate₀,growth_rate₁)
    Nₜ = @chain path begin
        readdlm(_) 
        get_cell_data(MachineState(),_) 
        get_population_competition_time_series(
                MachineState(),_,
                growth_rate₀,growth_rate₁)
    end
    return Nₜ
end