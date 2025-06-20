include("RedfieldThreeLevelEngine.jl")


savePlotsAt = "/home/gautham/ThreeLevelEngineWithQuantumLoad/FiguresSimulation/RedfieldCouplingStrength/"
saveTimeEvolDataAt = "/home/gautham/ThreeLevelEngineWithQuantumLoad/TimeEvolData/RedfieldCouplingStrength/"


############### Plotting the mean and the variance rates of the load as a function of g ###############
using JLD2, Plots

mean_arr = []
var_arr = []

βh, ωl = 0.40, 1.0
g_arr = range(0.0005, 0.75, 40)
δ_arr = range(-0.2, 0.2, 5)
δ_arr = [-0.5, -0.25, 0, 0.25, 0.5]

# Iterate over the coupling strength if the mean-variance data does not exist
for δ in δ_arr
    means_g = []
    vars_g = []

    for g in g_arr
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

        push!(means_g, mean)
        push!(vars_g, var)
        println("Loaded the data for β=$βh, ω=$ωl, g=$g, δ=$δ.")
    end

    # println("Means: ", means)
    push!(mean_arr, means_g)
    push!(var_arr, vars_g)
end

# Plot the mean and variance rates
MeanPlot = plot(
    xlabel=L"g", 
    ylabel=L"\dot{⟨n_l⟩}", 
    title="Mean rate for βh, ωl, δ= $(βh), $(ωl), $δ"
)
VarPlot = plot(
    xlabel=L"g", 
    ylabel=L"⟨Δ\dot{n}_l⟩^2", 
    title="Variance rate for βh, ωl, δ= $(βh), $(ωl), $δ"
)

# Generate labels for each plot and reshape them into row matrices
mean_labels = reshape(["δ = $δ" for δ in δ_arr], 1, :)
var_labels = reshape(["δ = $δ" for δ in δ_arr], 1, :)

plot!(
    MeanPlot, 
    g_arr, 
    mean_arr, 
    labels=mean_labels, 
    lw=2
)
plot!(
    VarPlot, 
    g_arr, 
    var_arr, 
    labels=var_labels, 
    lw=2
)

savefig(plot(MeanPlot, VarPlot, layout=(2, 1), size=(600, 800), margin=5mm), "RedfieldLoadMeanVarianceRate_β_$(βh)_ω_$(ωl).svg") |> println;
