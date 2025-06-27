# Use parallelization techniques to compute the dynamics info for all parameter values in our phase space
using Distributed
addprocs(40)

@everywhere include("RedfieldThreeLevelEngine.jl")
println("Loaded the RedfieldThreeLevelEngine module.")

savePlotsAt = "/home/gautham/ThreeLevelQuantumEngine/RedfieldDetuningNCoupling/Figures/"
saveTimeEvolDataAt = "/home/gautham/ThreeLevelQuantumEngine/RedfieldDetuningNCoupling/TimeEvolData/"


############### Plotting the mean and the variance rates of the load as a function of g ###############
using Plots, IterTools

mean_arr = []
var_arr = []

βh, ωl = 0.25, 3.0

parameter="δ" # "g" or "δ"

# Set the arrays depending on what parameter we are sweeping
if parameter == "g"
    # For the coupling strength
    g_arr = range(0.0, 0.4, 100)
    δ_arr = range(-0.5, 0.5, 5)
    δ_arr = [-0.5, -0.2, 0.001, 0.2, 0.5]

elseif parameter == "δ"
    # For the detuning for different coupling strengths
    g_arr = [0.0005, 0.005, 0.05, 0.15, 0.3]
    δ_arr = range(-1.5, 1.5, 100)
    # include 0 detuning in the plot
    δ_arr = [δ_arr[1:length(δ_arr)÷2]..., 0, δ_arr[length(δ_arr)÷2+1:end]...]
else
    println("Invalid parameter. Exiting...")
    exit()
end

# The folders where we save the data and the figures
if !isdir("TimeEvolData/$parameter")
    mkdir("TimeEvolData/$parameter")
end

if !isdir("Figures/$parameter")
    mkdir("Figures/$parameter")
end

saveTimeEvolDataAt = "TimeEvolData/$parameter/"
savePlotsAt = "Figures/$parameter/"

# Iterate over the coupling strength if the mean-variance data does not exist
@time @sync @distributed for (δ, g) in collect(IterTools.product(δ_arr, g_arr))
    # Set the filename to save the data
    filename = "$(saveTimeEvolDataAt)/RedfieldThreeLevelEngine_β_$(βh)_ω_$(ωl)_g_$(g)_δ_$δ.jld2"

    # Check if we already have the data
    if !isfile(filename)
        println("Data does not exist for β=$βh, ω=$ωl, g=$g, δ=$δ. Running the simulation...")
        local ergotropy_rate, energy_rate, means, variances, Qh, Qc, t = main(βh, ωl, g, δ; steadyStateDynamics=true, savePlots=true, savePlotsAt=savePlotsAt)
        println("Qh, Qc: ", Qh, ", ", Qc)
        @save filename ergotropy_rate energy_rate  means variances Qh Qc t
        println("Ran the simulation for β=$βh, ω=$ωl, g=$g, δ=$δ and saved the data.")
    end
end

println("Loading the data...")

for δ in δ_arr
    means_g = []
    vars_g = []

    for g in g_arr
        # Set the filename to save the data
        filename = "$(saveTimeEvolDataAt)RedfieldThreeLevelEngine_β_$(βh)_ω_$(ωl)_g_$(g)_δ_$δ.jld2"

        # Load the data from the jld2 file for the given βh and ωl
        println("Loading the data...")
        
        local ergotropy_rate, energy_rate, means, variances, t
        @load filename ergotropy_rate energy_rate means variances t
        local n = length(t)

        t, means, variances = t[n÷2:end], means[n÷2:end], variances[n÷2:end]
        println("Length of means: ", length(means), ", Length of variances: ", length(variances))
        local mean, var = linreg(t, means)[2], linreg(t, variances)[2]

        push!(means_g, mean)
        push!(vars_g, var)
    end

    # println("Means: ", means)
    push!(mean_arr, means_g)
    push!(var_arr, vars_g)
end

println("Loaded all the data. Plotting the phase spaces...")

# convert the mean and variance arrays to matrices
mean_arr = hcat(mean_arr...)
var_arr = hcat(var_arr...)

if parameter == "g"
    println("Plotting the n vs g for different δ...")
    # Plot the mean and variance rates
    MeanPlot = plot(
        xlabel=L"g", 
        ylabel=L"\dot{⟨n_l⟩}", 
        xlim=(minimum(g_arr), maximum(g_arr)),
        title="Mean rate for βh= $(βh), ωl $(ωl)"
    )
    VarPlot = plot(
        xlabel=L"g", 
        ylabel=L"⟨Δ\dot{n}_l⟩^2", 
        xlim=(minimum(g_arr), maximum(g_arr)),
        title="Variance rate for βh= $(βh), ωl $(ωl)"
    )

    # Generate labels for each plot and reshape them into row matrices
    δ_labels = reshape(["δ = $δ" for δ in δ_arr], 1, :)
    
    hline!(MeanPlot, [0], lw=2, label=L"\dot{⟨n_l⟩}=0", linestyle=:dash, color=:black)
    plot!(
        MeanPlot, 
        g_arr, 
        mean_arr, 
        labels=δ_labels, 
        lw=2
    )
    plot!(
        VarPlot, 
        g_arr, 
        var_arr, 
        labels=δ_labels, 
        lw=2
    )

elseif parameter == "δ"
    mean_arr = mean_arr'
    var_arr = var_arr'
    
    println("Plotting the n vs δ for different g...")
    # Plot the mean and variance rates
    MeanPlot = plot(
        xlabel=L"δ", 
        ylabel=L"\dot{⟨n_l⟩}", 
        xlim=(minimum(δ_arr), maximum(δ_arr)), 
        yformatter=(x -> x * 1e4),
        size=(400, 250), 
        top_margin=5mm,
        # title="Mean rate for βh= $(βh), ωl $(ωl)"
    )
    VarPlot = plot(
        xlabel=L"δ", 
        ylabel=L"⟨Δ\dot{n}_l⟩^2", 
        xlim=(minimum(δ_arr), maximum(δ_arr)), 
        yformatter = (x -> x * 1e4),
        # title="Variance rate for βh= $(βh), ωl $(ωl)"
    )

    # Generate labels for each plot and reshape them into row matrices
    g_labels = reshape(["g = $g" for g in g_arr], 1, :)

    # hline!(MeanPlot, [0], lw=2, label="", linestyle=:dash, color=:black)
    plot!(
        MeanPlot, 
        δ_arr, 
        mean_arr, 
        labels=g_labels, 
        legendfontsize=7,
        lw=2
    )
    annotate!(MeanPlot, -1.3, .00052, text(L"× 10^{-4}", :center, 12, "black"))
    plot!(
        VarPlot, 
        δ_arr, 
        var_arr, 
        labels=g_labels, 
        lw=2
    )
    annotate!(VarPlot, -1.5, .0031, text(L"× 10^{-4}", :center, 12, "black"))

end

savefig(plot(MeanPlot, VarPlot, layout=(1, 2), size=(850, 350), margin=5mm), "RedfieldLoadMeanVarianceRate_β_$(βh)_ω_$(ωl)_for_$parameter.svg") |> println;

savefig(MeanPlot, "RedfieldJustDetuning_β_$(βh)_ω_$(ωl)_for_$parameter.svg") |> println;
