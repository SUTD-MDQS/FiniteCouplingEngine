# use ThreeLevelEngineWithQuantumLoad.jl to run the time evolution and then use the biased diffusion in the load to plot
include("/home/gautham/ThreeLevelEngineWithQuantumLoad/LocalLindbladDynamics/ThreeLevelEngineWithQuantumLoad.jl")
using Plots.PlotMeasures

# the system parameters
# the sweep parameters that we vary to enter different regimes
βh, ωl = 0.5, 1.0
g = 0.05
tf = 2000

ω = [1, ωl]
β = [βh, 2.0]
el, l = [2, 40]

# define the fock basis
engineBasis = FockBasis(el)
loadBasis = FockBasis(l)

# run a simulation while returning the solution
_, _, _, means, variances, t, sol = main(βh, ωl, g; tf=tf, returnSol=true, steadyStateDynamics=true, savePlots=false)

# calculate the reduced density matrix
_, load_populations = calculate_populations(sol, engineBasis, loadBasis)


color_scale = cgrad([:white, :blue, :red])
mean_col = "#1B5E20"
var_col = "#D55E00"

# The biased diffusion in the load
Fig1a = heatmap(t, 0:loadBasis.N, load_populations, 
    xlabel=L"t", ylabel=L"p(n)", 
    size=(600, 300),
    color=color_scale, 
    colorbar=:none,
    legend=false,
    xlim=(0, t[end]),
    ylim=(12, loadBasis.N),
    xtickfontsize=20,
    ytickfontsize=20,
    guidefontsize=25,
    framestyle=:box,
)

t_intervals = 35
# Plot the variance bars at specific time intervals
t_bar = t[1:t_intervals:end]  # Time intervals
mean_bar = (means[1:t_intervals:end] .- means[1] .+ 20)  # Mean values
error_bar = sqrt.(variances[1:t_intervals:end])  # Error values (define this array)

# Plot the mean values
plot!(Fig1a, 
    t, (means .- means[1] .+ 20), 
    color=mean_col, 
    linestyle=:dashdot, 
    lw=3,
)

# Add scatter plot with error bars
plot!(Fig1a, 
    t_bar, mean_bar, 
    seriestype=:scatter, 
    markersize=5, 
    color=mean_col, 
    yerror=error_bar, 
    markerstrokewidth=3, 
    markerstrokecolor=var_col,
)

pyplot()
Fig1b = plot(
    t, means, 
    color=mean_col, 
    linewidth=3, 
    label=L"\mu(t)", 
    xlabel=L"t", 
    ylabel="", 
    xlim=(0, t[end]),
    xtickfontsize=20,
    ytickfontsize=20,
    guidefontsize=25,
    legendfontsize=25,
    size=(500,500), 
)

# Add the second plot (variance with right y-axis)
plot!(
    twinx(),
    t, variances, 
    color=var_col, 
    linewidth=3, 
    label=L"\sigma(t)", 
    ylabel="", 
    xlim=(0, t[end]),
    ylim=1.5 .* (minimum(variances),maximum(variances)),
    xticks=:none, 
    ytickfontsize=20,
    guidefontsize=25,
    legendfontsize=25,
    framestyle=:box,
)
layout = @layout [a{0.65w} b{0.35w}]
Fig1 = plot(Fig1a, Fig1b, layout=layout, size=(1700, 600), margins=12mm)

savefig(Fig1, "BiasedDiffusionWithMeanAndVar.svg")
