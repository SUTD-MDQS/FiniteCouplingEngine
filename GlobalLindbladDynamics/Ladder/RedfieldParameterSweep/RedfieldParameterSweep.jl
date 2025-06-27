# Use parallelization techniques to compute the dynamics info for all parameter values in our phase space
using Distributed
addprocs(40)

@everywhere include("/home/gautham/ThreeLevelLadderEngine/RedfieldParameterSweep/RedfieldThreeLevelLadder.jl")
println("Loaded the RedfieldThreeLevelLadder module.")

savePlotsAt = "/home/gautham/ThreeLevelLadderEngine/RedfieldParameterSweep/Figures/"
saveTimeEvolDataAt = "/home/gautham/ThreeLevelLadderEngine/RedfieldParameterSweep/TimeEvolData/"


############### Plotting the mean and the variance rates of the load as a function of g ###############
using Plots, IterTools

mean_arr = []
var_arr = []
regimes_arr = []

βh_arr = range(0.01, 2.0, 60)
δ, ωl = 0.0, 3.0
g_arr = range(0.0, 1.5, 60)

# δ_arr = range(-0.2, 0.2, 5)
# δ_arr = [-0.5, -0.25, 0, 0.25, 0.5]


# Parallelize the data computation and saving
@time @sync @distributed for (βh, g) in collect(IterTools.product(βh_arr, g_arr))
    # Set the filename to save the data
    filename = "$(saveTimeEvolDataAt)RedfieldThreeLevelEngine_β_$(βh)_ω_$(ωl)_g_$(g)_δ_$δ.jld2"

    # Check if we already have the data
    if !isfile(filename)
        println("Data does not exist for β=$βh, ω=$ωl, g=$g, δ=$δ. Running the simulation...")
        local ergotropy_rate, energy_rate, means, variances, Qh, Qc, t = main(βh, ωl, g, δ; steadyStateDynamics=true, savePlots=true, savePlotsAt=savePlotsAt)

        @save filename ergotropy_rate energy_rate means variances Qh Qc t
        println("Ran the simulation for β=$βh, ω=$ωl, g=$g, δ=$δ and saved the data.")
    end
end

println("Found all the data. Loading the data...")

# Serial loop for loading and processing the data
for βh in βh_arr
    means_g = []
    vars_g = []
    regimes_g = []

    for g in g_arr
        # Set the filename to load the data
        filename = "$(saveTimeEvolDataAt)RedfieldThreeLevelEngine_β_$(βh)_ω_$(ωl)_g_$(g)_δ_$δ.jld2"

        # Load the data from the jld2 file for the given βh and ωl
        println("Loading the data for β=$βh, ω=$ωl, g=$g, δ=$δ...")
        local ergotropy_rate, energy_rate, means, variances, Qh, Qc, t

        @load filename ergotropy_rate energy_rate means variances Qh Qc t
        local n = length(t)
        
        t, means, variances = t[n÷2:end], means[n÷2:end], variances[n÷2:end]
        println("Length of means: ", length(means), ", Length of variances: ", length(variances))
        local mean, var = linreg(t, means)[2], linreg(t, variances)[2]

        push!(means_g, mean)
        push!(vars_g, var)

        # Determine the regime of the system
        if energy_rate >= 0 && Qh > 0 && Qc < 0
            push!(regimes_g, "engine")
        elseif energy_rate < 0 && Qh >= 0 && Qc < 0
            push!(regimes_g, "accelerator")
        elseif energy_rate < 0 && Qh < 0 && Qc >= 0
            push!(regimes_g, "refrigerator")
        elseif energy_rate < 0 && Qh < 0 && Qc < 0
            push!(regimes_g, "heater")
        else
            energy_exp = energy_rate == 0 ? 0 : Int(floor(log10(abs(energy_rate))))
            Qh_exp = Qh == 0 ? 0 : Int(floor(log10(abs(Qh))))
            Qc_exp = Qc == 0 ? 0 : Int(floor(log10(abs(Qc))))
            push!(regimes_g, "unknown. E:$(Int(sign(energy_rate)))e$energy_exp, Qh:$(Int(sign(Qh)))e$Qh_exp, Qc:$(Int(sign(Qc)))e$Qc_exp")
        end

        println("Loaded the data for β=$βh, ω=$ωl, g=$g, δ=$δ.")
    end

    push!(mean_arr, means_g)
    push!(var_arr, vars_g)
    push!(regimes_arr, regimes_g)
end

println("Loaded all the data. Plotting the phase spaces...")

# convert the mean and variance arrays to matrices
mean_arr = hcat(mean_arr...)'
var_arr = hcat(var_arr...)'

min_mean, max_mean = minimum(mean_arr), maximum(mean_arr)
min_var, max_var = minimum(var_arr), maximum(var_arr)

zero_pos_mean = (0 - min_mean) / (max_mean - min_mean)
zero_pos_var = (0 - min_var) / (max_var - min_var)

rdbu_mean = cgrad([:blue, :white, :red], [0, zero_pos_mean, 1])
rdbu_var = cgrad([:white, :blue, :red], [0, zero_pos_var, 1])


# convert the regimes array to a matrix of contours for the heatmap. engine = 1, accelerator = 2, refrigerator = 3, heater = 4, unknown = 0
regimes_arr = hcat(regimes_arr...)
regimes_arr = map(x -> x == "engine" ? 1 : x == "accelerator" ? 2 : x == "refrigerator" ? 3 : x == "heater" ? 4 : 0, regimes_arr)'

engine_matrix = regimes_arr .== 1 |> Int
accelerator_matrix = regimes_arr .== 2 |> Int
refrigerator_matrix = regimes_arr .== 3 |> Int
heater_matrix = regimes_arr .== 4 |> Int

engine_matrix *= max_mean
accelerator_matrix *= max_mean
refrigerator_matrix *= max_mean
heater_matrix *= max_mean

# Plot the mean and variance rates
MeanPlot = heatmap(
    g_arr[1:end], 
    βh_arr[1:end],
    mean_arr[1:end,1:end], 
    xlabel=L"g", 
    ylabel=L"\beta_h", 
    xlims=(g_arr[1], g_arr[end]),
    ylims=(βh_arr[1], βh_arr[end]),
    # L"\dot{⟨n_l⟩}"
    title=L"\dot{⟨n_l⟩}"*" for ωl, δ= $(ωl), $δ", 
    color=rdbu_mean
)
# VarPlot = heatmap(
#     g_arr[2:end], 
#     βh_arr[2:end],
#     var_arr[2:end,2:end],
#     xlabel=L"g", 
#     ylabel=L"\beta_h", 
#     # L"⟨Δ\dot{n}_l⟩^2"
#     title=L"⟨Δ\dot{n}_l⟩^2"*" for ωl, δ= $(ωl), $δ", 
#     color=rdbu_var
# )
contour!(
    MeanPlot,
    g_arr, 
    βh_arr,
    mean_arr, 
    levels=[0]
)
# contour!(MeanPlot, g_arr, βh_arr, engine_matrix, levels=[0.00005], linewidth=2, color=:black, label="Engine")
contour!(MeanPlot, g_arr, βh_arr, accelerator_matrix, levels=[0.00005], linewidth=2, color=:red, label="Accelerator")
contour!(MeanPlot, g_arr, βh_arr, refrigerator_matrix, levels=[0.00005], linewidth=2, color=:blue, label="Refrigerator")
contour!(MeanPlot, g_arr, βh_arr, heater_matrix, levels=[0.00005], linewidth=2, color=:green, label="Heater")
println("Saving the plot...")
savefig(plot(MeanPlot, layout=(1, 1), size=(1000, 1000), margin=5mm), "RedfieldLoadMeanRate_ω_$(ωl)_δ_$δ.svg") |> println;

# using GLMakie
# fig = Figure(resolution = (600, 1000))

# # Plot the mean rates heatmap
# ax1 = Axis(fig[1, 1], title = L"\dot{⟨n_l⟩}"*" for ωl, δ= $(ωl), $δ", xlabel = L"g", ylabel = L"\beta_h")
# heatmap!(ax1, g_arr, βh_arr, mean_arr, colormap = rdbu_mean)
# contour!(ax1, g_arr, βh_arr, mean_arr, levels = [0], color = :black)
# contour!(ax1, g_arr, βh_arr, accelerator_matrix, levels = [0.00005], linewidth = 2, color = :red, label = "Accelerator")
# contour!(ax1, g_arr, βh_arr, refrigerator_matrix, levels = [0.00005], linewidth = 2, color = :blue, label = "Refrigerator")
# contour!(ax1, g_arr, βh_arr, heater_matrix, levels = [0.00005], linewidth = 2, color = :green, label = "Heater")

# # Plot the variance rates heatmap
# ax2 = Axis(fig[2, 1], title = L"⟨Δ\dot{n}_l⟩^2"*" for ωl, δ= $(ωl), $δ", xlabel = L"g", ylabel = L"\beta_h")
# heatmap!(ax2, g_arr, βh_arr, var_arr, colormap = rdbu_var)

# # Save the plots
# println("Saving the plot...")
# save("MeanPlot.png", fig)
# save("VarPlot.png", fig)
# save("RedfieldLoadMeanVarianceRate_ω_$(ωl)_δ_$δ.png", fig)
# println("Plot saved.")
