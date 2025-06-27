include("/home/gautham/ThreeLevelEngineWithQuantumLoad/LocalLindbladDynamics/ThreeLevelEngineWithQuantumLoad.jl")

savePlotsAt = "/home/gautham/ThreeLevelEngineWithQuantumLoad/FiguresSimulation/CouplingStrength/"
saveTimeEvolDataAt = "/home/gautham/ThreeLevelEngineWithQuantumLoad/TimeEvolData/CouplingStrength/"


############### Plotting the mean and the variance rates of the load as a function of g ###############
using JLD2, Plots

mean_arr = []
var_arr = []

δ, ωl = 0.0, 1.0
β_arr = range(0.5, 1.75, 6)
β_arr = [0.5, 0.75, 1.0, 1.25, 1.75]
g_arr = range(0.0005, 0.1, 50)
# g_arr = vcat(g_left, g_right)

# Iterate over the coupling strength if the mean-variance data does not exist
for βh in β_arr
    means_g = []
    vars_g = []

    for g in g_arr
        # Set the filename to save the data
        filename = "$(saveTimeEvolDataAt)ThreeLevelEngine_β_$(βh)_ω_$(ωl)_g_$(g)_δ_$δ.jld2"

        # Check if we already have the data
        if !isfile(filename)
            println("Data does not exist for β=$βh, ω=$ωl, g=$g, δ=$δ. Running the simulation...")
            local ergotropy_rate, energy_rate, pop_inversion, means, variances, t, Qh, Qc = main(βh, ωl, g, δ; steadyStateDynamics=true, savePlots=true)

            @save filename ergotropy_rate energy_rate pop_inversion means variances t Qh Qc
            println("Ran the simulation for β=$βh, ω=$ωl, g=$g, δ=$δ and saved the data.")
        end

        # Load the data from the jld2 file for the given βh and ωl
        println("Loading the data...")
        
        @load filename means variances t
        local mean, var = linreg(t, means)[2], linreg(t, variances)[2]

        push!(means_g, mean)
        push!(vars_g, var)
        println("Loaded the data for β=$βh, ω=$ωl, g=$g, δ=$δ.")
    end

    # println("Means: ", means)
    push!(mean_arr, means_g)
    push!(var_arr, vars_g)
end

# Generate labels for each plot and reshape them into row matrices
mean_labels = reshape([L"β_h = "*"$βh"*L"(\hbar\omega_c)^{-1}" for βh in β_arr], 1, :)

heat_clr = "#8c5b89"
eng_clr = "#3a873c"
acc_clr = "#e5b12e"
ref_clr = "#50b3d1"
unknown_clr = "#cccccc"

# Map regimes to colors
regime_colors = Dict(
    "engine" => eng_clr,
    "accelerator" => acc_clr,
    "heater" => heat_clr,
    "refrigerator" => ref_clr,
    "unknown" => unknown_clr
)
regime_markers = Dict(
    "engine" => :square,
    "accelerator" => :diamond,
    "heater" => :circle,
    "refrigerator" => :utriangle
)
line_colors = reshape(
    [
    "#ff99aaff",
    "#e50000ff", 
    "#631200ff", 
    "#4400dbff", 
    "#a1b8ffff"
    ], 
1,:
)


# Plot the mean and variance rates
MeanPlot = plot(
    xlabel=L"g", 
    ylabel=L"v", 
    xlim=(minimum(g_arr), maximum(g_arr)),
    # title="Mean rate for βh, ωl, δ= $(βh), $(ωl), $δ", 
    yformatter = (x -> x * 1e4),
    size=(500,400), 
    legendfontsize=7, 
    legend=:outerbottom,
    legendcolumns=3, 
    framestyle=:box,
)
plot!(
    MeanPlot, 
    g_arr, 
    mean_arr, 
    labels=mean_labels, 
    color=line_colors,
    lw=3, 
)
hline!(MeanPlot, 
    [0], 
    label=L"v=0", 
    linestyle=:dash, 
    linecolor=:black,
    lw=3, 
)
# Define the interval for selecting points
interval = 4  # Adjust this to control the number of points plotted

# Loop through the mean_arr and g_arr
for mean_values in mean_arr
    # Select points at equal intervals
    selected_indices = 1:interval:length(g_arr)
    selected_g = g_arr[selected_indices]
    selected_means = mean_values[selected_indices]

    # Classify points into regimes
    regimes = [v>0 ? "engine" : "refrigerator" for v in selected_means]

    # Map regimes to colors and markers
    scatter_colors = [regime_colors[regime] for regime in regimes]
    scatter_markers = [regime_markers[regime] for regime in regimes]

    # Plot the scatter points with regime-specific colors and markers
    scatter!(
        MeanPlot,
        selected_g, selected_means,
        color=scatter_colors,
        marker=scatter_markers,
        label="",
        markersize=5
    )
end

VarPlot = plot(
    xlabel=L"g", 
    ylabel=L"D", 
    xlim=(minimum(g_arr), maximum(g_arr)),
    # title="Variance rate for βh, ωl, δ= $(βh), $(ωl), $δ", 
    size=(500,400), 
    legendfontsize=7, 
    framestyle=:box,
)
plot!(
    VarPlot, 
    g_arr, 
    var_arr, 
    labels=mean_labels, 
    lw=2, 
)

# savefig(plot(MeanPlot, VarPlot, layout=(2, 1), size=(600, 800), margin=5mm), "MeanVarianceRate_β_$(βh)_ω_$(ωl).svg") |> println;

savefig(MeanPlot, "MeanVsG.svg")
# savefig(VarPlot, "VarVsG.svg")
