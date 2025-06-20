# this script runs the engine efficiency simulation and plots the results for a given βh across ωl
using JLD2
include("/home/gautham/ThreeLevelEngineWithQuantumLoad/LocalLindbladDynamics/ThreeLevelEngineWithQuantumLoad.jl")

# the parameters
βh = 1.0
g = 0.05
tf = 3000

# the range of ωl values
ωl_values = 0.1:0.2:5.5

ηc = 1 - 1/4

# the simulation
for i in eachindex(ωl_values)
    ωl = ωl_values[i]
    ergotropy_rate, energy_rate, pop_inversion, means, variances, t, Qh, Qc = main(βh, ωl, g; tf=tf, steadyStateDynamics=true)
    # if the engine draws energy from the hot bath, calculate the efficiency else set it to zero
    Qh > Qc ? η = energy_rate / Qh : η = 0
    println("η = $η, ηc = $ηc")

    # save the results
    save("TimeEvolData/EngineEfficiency/efficiency_ω_$(ωl)_ηc_$(ηc).jld2","ηc", ηc, "η", η, "ωl", ωl)
end

# load the data into an array
η_values = zeros(length(ωl_values))

for i in 1:length(ωl_values)
    data = load("TimeEvolData/EngineEfficiency/efficiency_ω_$(ωl_values[i])_ηc_$(ηc).jld2")
    η_values[i] = data["η"]
end

# plot the results
Fig = plot(
    xlabel=L"ω_e", 
    ylabel=L"η_e", 
    xlim=(ωl_values[1], ωl_values[end]),
    ylim=(0, 1),
    size=(400, 400), 
    framestyle=:box,
    legend=:bottomright,
)
plot!(Fig, 
    ωl_values, 
    η_values, 
    label=L"η_e", 
    linestyle=:solid, 
    marker=:circle, 
    markersize=2, 
)
hline!([ηc], label=L"η_c", linestyle=:dash)
plot!(Fig, 
    ωl_values, 
    ωl_values ./ (ωl_values .+ 1),
    linestyle=:dashdot, 
    label=L"η = \frac{ω_e}{ω_e + 1}",
)
