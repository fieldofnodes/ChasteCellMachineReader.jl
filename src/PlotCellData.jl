"""
    For each column in the cell data frame,
    we plot the corresponding time series as well the per capita version
"""
function plot_cell_populations_TS(cdf::DataFrame)
    tim = cdf.time
    att = cdf.attacker
    tar = cdf.target
    m₁ = cdf.m₁
    m₂ = cdf.m₂
    m₃ = cdf.m₃
    mₜ = cdf.mₜ
    m₁_att = m₁ ./ att
    m₂_att = m₂ ./ att
    m₃_att = m₃ ./ att
    mₜ_att = mₜ ./ att
    n_d_att = cdf.neigh_diff_to_attacker
    n_d_tar = cdf.neigh_diff_to_target
    n_d_att_p_cap = n_d_att ./ att
    n_d_tar_p_cap = n_d_tar ./ tar


    # Plot target and attackers
    f₁ = Figure(resolution=(800,600))
    ax = Axis(f₁[1,1],xlabel="Time",ylabel="Count")
    lines!(ax,tim,tar,label="Target")
    lines!(ax,tim,att,label="Attacker")
    axislegend("Cell type")
    f₁

    # Plot machines
    f₂ = Figure(resolution=(800,600))
    ax = Axis(f₂[1,1],xlabel="Time",ylabel="Count")
    lines!(ax,tim,m₁,label="State 1")
    lines!(ax,tim,m₂,label="State 2")
    lines!(ax,tim,m₃,label="State 3")
    axislegend("Machine states")
    f₂


    # Plot per capita machines
    f₃ = Figure(resolution=(800,600))
    ax = Axis(f₃[1,1],xlabel="Time",ylabel="Per capita count")
    lines!(ax,tim,m₁_att,label="State 1")
    lines!(ax,tim,m₂_att,label="State 2")
    lines!(ax,tim,m₃_att,label="State 3")
    axislegend("Per capita \nMachine states")
    f₃


    # Plot total machines
    f₄ = Figure(resolution=(800,600))
    ax = Axis(f₄[1,1],xlabel="Time",ylabel="Count")
    lines!(ax,tim,mₜ,label="Total")
    axislegend("Machines")
    f₄

    # Plot total per capita machines
    f₅ = Figure(resolution=(800,600))
    ax = Axis(f₅[1,1],xlabel="Time",ylabel="Per capita count")
    lines!(ax,tim,mₜ_att,label="Total Machines")
    axislegend("Per capita \nMachines")
    f₅


    # Plot total nerighbours
    f₆ = Figure(resolution=(800,600))
    ax = Axis(f₆[1,1],xlabel="Time",ylabel="Count")
    lines!(ax,tim,n_d_att,label="Attacker")
    lines!(ax,tim,n_d_tar,label="Target")
    axislegend("Neighbours \ndifferent to")
    f₆

    # Plot total per capita neighbours
    f₇ = Figure(resolution=(800,600))
    ax = Axis(f₇[1,1],xlabel="Time",ylabel="Per capita count")
    lines!(ax,tim,n_d_att_p_cap,label="Attacker")
    lines!(ax,tim,n_d_tar_p_cap,label="Target")
    axislegend("Per capita \nNeighbours \ndifferent to")
    f₇

    return [f₁,f₂,f₃,f₄,f₅,f₆,f₇]
end
