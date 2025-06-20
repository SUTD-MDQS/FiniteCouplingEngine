include("RedfieldThreeLevelEngine.jl")

savePlotsAt = "/home/gautham/ThreeLevelEngineWithQuantumLoad/FiguresSimulation/RedfieldDetuning/"
saveTimeEvolDataAt = "/home/gautham/ThreeLevelEngineWithQuantumLoad/TimeEvolData/RedfieldDetuning/"
############### Plotting the mean and the variance rates of the load as a function of g ###############
using JLD2, Plots

mean_arr = []
var_arr = []

βh, ωl = 0.40, 1.0
δ_arr = range(-5, 5, 50)
g_arr = [0.0005, 0.005, 0.05, 0.1, 0.75]

# Iterate over the coupling strength if the mean-variance data does not exist
for g in g_arr
    means_δ = []
    vars_δ = []

    for δ in δ_arr
        # Set the filename to save the data
        filename = "$(saveTimeEvolDataAt)RedfieldThreeLevelEngine_β_$(βh)_ω_$(ωl)_g_$(g)_δ_$δ.jld2"

        # Check if we already have the data
        if !isfile(filename)
            println("Data does not exist for β=$βh, ω=$ωl, g=$g, δ=$δ. Running the simulation...")
            local ergotropy_rate, energy_rate, means, variances, t = main(βh, ωl, g, δ; steadyStateDynamics=true, savePlots=true)

            @save filename ergotropy_rate energy_rate  means variances t
            println("Ran the simulation for β=$βh, ω=$ωl, g=$g, δ=$δ and saved the data.")
        end

        # Load the data from the jld2 file for the given βh and ωl
        println("Loading the data...")
        
        local ergotropy_rate, energy_rate, means, variances, t
        @load filename ergotropy_rate energy_rate means variances t

        local mean, var = linreg(t, means)[2], linreg(t, variances)[2]

        push!(means_δ, mean)
        push!(vars_δ, var)
        println("Loaded the data for β=$βh, ω=$ωl, g=$g, δ=$δ.")
    end

    # println("Means: ", means)
    push!(mean_arr, means_δ)
    push!(var_arr, vars_δ)
end


# Plot the mean and variance rates
MeanPlot = plot(
    xlabel=L"\delta", 
    ylabel=L"\dot{⟨n_l⟩}", 
    title="Mean rate for βh, ωl, g= $(βh), $(ωl), $g"
)
VarPlot = plot(
    xlabel=L"\delta", 
    ylabel=L"⟨Δ\dot{n}_l⟩^2", 
    title="Variance rate for βh, ωl, g= $(βh), $(ωl), $g"
)

# Generate labels for each plot and reshape them into row matrices
mean_labels = reshape(["g = $g" for g in g_arr], 1, :)
var_labels = reshape(["g = $g" for g in g_arr], 1, :)

plot!(
    MeanPlot, 
    δ_arr, 
    mean_arr, 
    labels=mean_labels, 
    lw=2
)
plot!(
    VarPlot, 
    δ_arr, 
    var_arr, 
    labels=var_labels, 
    lw=2
)

savefig(plot(MeanPlot, VarPlot, layout=(2, 1), size=(600, 800), margin=5mm), "RedfieldDetuningRate_β_$(βh)_ω_$(ωl)_g_$g.svg") |> println;
