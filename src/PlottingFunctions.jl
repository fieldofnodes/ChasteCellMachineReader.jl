
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
