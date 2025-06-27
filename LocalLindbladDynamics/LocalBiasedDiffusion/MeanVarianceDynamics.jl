include("ThreeLevelEngineWithQuantumLoad.jl")

############# Plotting the mean and the variance of the load #############
using JLD2, Plots

FigMean = plot(xlabel="Time", ylabel="Mean", title="Mean of Load Populations for βh, ωl = $(βh), $(ωl)");
FigVar = plot(xlabel="Time", ylabel="Variance", title="Variance of Load Populations for βh, ωl = $(βh), $(ωl)");

mean_arr = []
var_arr = []

βh, ωl = 0.45, 1.0
g_arr = range(0.0005, 0.1, 70)

# Iterate over the coupling strength if the mean-variance data does not exist
for g in g_arr
    filename = "$(savePlotsAt)TimeEvolData/CouplingStrength/ThreeLevelEngineWithQuantumLoad_β_$(βh)_ω_$(ωl)_g_$(g).jld2"
    if !isfile(filename)
        ergotropy_rate, energy_rate, pop_inversion, means, variances, t = main(βh, ωl, g; steadyStateDynamics=true)

        plot!(FigMean, t, means, label="g = $g")
        plot!(FigVar, t, variances, label="g = $g")
        mean = linreg(t, means)[2]
        push!(mean_arr, mean)
        var = linreg(t, variances)[2]
        push!(var_arr, var)
        
        # save the values in a jld2 file
        @save filename ergotropy_rate energy_rate pop_inversion means variances t
    end
end

# Load the data from the jld2 file for the given βh and ωl
for g in g_arr
    filename = "$(savePlotsAt)TimeEvolData/CouplingStrength/ThreeLevelEngineWithQuantumLoad_β_$(βh)_ω_$(ωl)_g_$(g).jld2"
    @load filename ergotropy_rate energy_rate pop_inversion means variances t
    push!(mean_arr, linreg(t, means)[2])
    push!(var_arr, linreg(t, variances)[2])
end

# savefig(plot(FigMean, FigVar, layout=(2, 1), size=(600, 800), margin=2mm), "$(savePlotsAt)/LoadMeanVariance_β_$(βh)_ω_$(ωl).svg")

MeanPlot = plot(g_arr, mean_arr, xlabel=L"g", ylabel=L"\dot{⟨n_l⟩}", title="Mean rate for βh, ωl = $(βh), $(ωl)", label="Numerical value");
VarPlot = plot(g_arr, var_arr, xlabel=L"g", ylabel=L"⟨Δ\dot{n}_l⟩^2", title="Variance rate for βh, ωl = $(βh), $(ωl)", label="Numerical value");

savefig(plot(MeanPlot, VarPlot, layout=(2, 1), size=(600, 800), margin=5mm), "$(savePlotsAt)/LoadMeanVarianceRate_β_$(βh)_ω_$(ωl).svg") |> println;
