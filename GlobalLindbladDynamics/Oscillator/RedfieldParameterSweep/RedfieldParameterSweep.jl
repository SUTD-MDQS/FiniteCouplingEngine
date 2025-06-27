# Use parallelization techniques to compute the dynamics info for all parameter values in our phase space
using Distributed
addprocs(40)

@everywhere include("RedfieldThreeLevelEngine.jl")
println("Loaded the RedfieldThreeLevelEngine module.")

savePlotsAt = "/home/gautham/ThreeLevelQuantumEngine/RedfieldParameterSweep/Figures/"
saveTimeEvolDataAt = "/home/gautham/ThreeLevelQuantumEngine/RedfieldParameterSweep/TimeEvolData/"


############### Plotting the mean and the variance rates of the load as a function of g ###############
using Plots, IterTools, Contour

# The arrays of data we will save for the parameter sweep
mean_arr = []
var_arr = []
regimes_arr = []
Pl_arr = []
Qh_arr = []
Qc_arr = []

# Define the parameters for the sweep
βh_arr = range(0.01, 2.0, 60)
δ, ωe = 0.0, 3.0
g_arr = range(0.0, 0.4, 60)

# δ_arr = range(-0.2, 0.2, 5)
# δ_arr = [-0.5, -0.25, 0, 0.25, 0.5]


# Parallelize the data computation and saving
@time @sync @distributed for (βh, g) in collect(IterTools.product(βh_arr, g_arr))
    # Set the filename to save the data
    filename = "$(saveTimeEvolDataAt)RedfieldThreeLevelEngine_β_$(βh)_ω_$(ωe)_g_$(g)_δ_$δ.jld2"

    # Check if we already have the data
    if !isfile(filename)
        println("Data does not exist for β=$βh, ω=$ωe, g=$g, δ=$δ. Running the simulation...")
        local ergotropy_rate, energy_rate, means, variances, Qh, Qc, t = main(βh, ωe, g, δ; steadyStateDynamics=true, savePlots=true, savePlotsAt=savePlotsAt)

        @save filename ergotropy_rate energy_rate means variances Qh Qc t
        println("Ran the simulation for β=$βh, ω=$ωe, g=$g, δ=$δ and saved the data.")
    end
end

println("Found all the data. Loading the data...")

# Serial loop for loading and processing the data
for βh in βh_arr
    # the arrays across the g values
    means_g = []
    vars_g = []
    regimes_g = []
    Pl_g = []
    Qh_g = []
    Qc_g = []

    for g in g_arr
        # Set the filename to load the data
        filename = "$(saveTimeEvolDataAt)RedfieldThreeLevelEngine_β_$(βh)_ω_$(ωe)_g_$(g)_δ_$δ.jld2"

        # Load the data from the jld2 file for the given βh and ωl
        println("Loading the data for β=$βh, ω=$ωe, g=$g, δ=$δ...")
        local ergotropy_rate, energy_rate, means, variances, Qh, Qc, t

        @load filename ergotropy_rate energy_rate means variances Qh Qc t
        local n = length(t)
        
        t, means, variances = t[n÷2:end], means[n÷2:end], variances[n÷2:end]
        println("Length of means: ", length(means), ", Length of variances: ", length(variances))
        local mean, var = linreg(t, means)[2], linreg(t, variances)[2]

        push!(means_g, mean)
        push!(vars_g, var)
        push!(Pl_g, ωe * mean)
        push!(Qh_g, Qh)
        push!(Qc_g, Qc)

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

        println("Loaded the data for β=$βh, ω=$ωe, g=$g, δ=$δ.")
    end

    push!(mean_arr, means_g)
    push!(var_arr, vars_g)
    push!(regimes_arr, regimes_g)
    push!(Pl_arr, Pl_g)
    push!(Qh_arr, Qh_g)
    push!(Qc_arr, Qc_g)
end

println("Loaded all the data. Plotting the phase spaces...")

# Function to draw lines along the x/y axis at val with spacing
function draw_lines!(plot_obj, xlims, ylims, orientation::Tuple{Symbol, Float64}, n::Int; color::Symbol=:black, markersize=2)
    # Extract orientation and value
    direction, val = orientation

    if direction == :y
        # Generate n points for a horizontal line at y = val
        x_points = range(xlims[1], xlims[2], length=n)
        y_points = fill(val, n)  # y is constant
        scatter!(plot_obj, x_points, y_points, color=color, marker=:circle, markersize=markersize, markerstrokewidth=0, label="")
    elseif direction == :x
        # Generate n points for a vertical line at x = val
        y_points = range(ylims[1], ylims[2], length=n)
        x_points = fill(val, n)  # x is constant
        scatter!(plot_obj, x_points, y_points, color=color, marker=:circle, markersize=markersize, markerstrokewidth=0, label="")
    else
        error("Invalid orientation. Use :x for horizontal lines or :y for vertical lines.")
    end
end

# convert the mean and variance arrays to matrices
mean_arr = hcat(mean_arr...)
var_arr = hcat(var_arr...)
Qh_arr = hcat(Qh_arr...)
Qc_arr = hcat(Qc_arr...)
Pl_arr = hcat(Pl_arr...)

# truncate the matrices to 2:end
mean_arr = mean_arr[2:end, 2:end]
var_arr = var_arr[2:end, 2:end]
Qc_arr = Qc_arr[2:end, 2:end]
Qh_arr = Qh_arr[2:end, 2:end]
g_arr = g_arr[2:end]
βh_arr = βh_arr[2:end]
Pl_arr = Pl_arr[2:end, 2:end]

# convert the regimes array to a matrix of contours for the heatmap. engine = 1, accelerator = 2, refrigerator = 3, heater = 4, unknown = 0
regimes_arr = hcat(regimes_arr...)
regimes_arr = map(x -> x == "engine" ? 1 : x == "accelerator" ? 2 : x == "refrigerator" ? 3 : x == "heater" ? 4 : 0, regimes_arr)'

engine_matrix = regimes_arr .== 1 |> Int
accelerator_matrix = regimes_arr .== 2 |> Int
refrigerator_matrix = regimes_arr .== 3 |> Int
heater_matrix = regimes_arr .== 4 |> Int

# normalise the βh_arr for plotting
βh_arr = βh_arr ./ maximum(βh_arr)

# Create the main phase plot
MeanPlot = plot(
    xlim=(g_arr[begin], g_arr[end]), 
    ylim=(βh_arr[begin], βh_arr[end]), 
    xlabel=L"g/\hbar\omega_c", 
    ylabel=L"\beta_h/\beta_c", 
    legend=false, 
    size=(400, 400),
    right_margin=3mm,
    xticks=range(0, 0.4, length=5),
    yticks=range(0, 1.0, length=5),
    xtickfontsize=16,
    ytickfontsize=16,
    guidefontsize=18,
    legendfontsize=10,
)
# The color palette of our phase plot
heat_clr = "#8c5b89"
eng_clr = "#3a873c"
acc_clr = "#e5b12e"
ref_clr = "#50b3d1"

# 1. Cover the entire plot area with heater color
hspan!([βh_arr[1], βh_arr[end]], color=heat_clr, alpha=1)

# 2. Contour plot of Refrigerator
contour_result_Qc = Contour.contours(g_arr, βh_arr, Qc_arr, [0])  # Compute contours for Qc
paths_Qc = [Contour.coordinates(curve) for curve in Contour.lines(Contour.levels(contour_result_Qc)[1])]
for (x_coords, y_coords) in paths_Qc
    # Plot the contour line
    plot!(x_coords, y_coords, color=:black, lw=2, label="")

    # Fill the negative area below the contour line
    plot!(x_coords, y_coords, fillrange=maximum(βh_arr), color=ref_clr, alpha=1.0, lw=0, label="")
end

# 3. Contour plot of Accelerator
contour_result_Qh = Contour.contours(g_arr, βh_arr, Qh_arr, [0])  # Compute contours for Qh
paths_Qh = [Contour.coordinates(curve) for curve in Contour.lines(Contour.levels(contour_result_Qh)[1])]
for (x_coords, y_coords) in paths_Qh
    # Plot the contour line
    plot!(x_coords, y_coords, color=:black, lw=2, label="")

    # Fill the negative area below the contour line
    plot!(x_coords, y_coords, fillrange=minimum(βh_arr), color=acc_clr, alpha=1.0, lw=0, label="")
end

# 4. Contour plot of Engine
contour_result_mean = Contour.contours(g_arr, βh_arr, Pl_arr, [0])  # Compute contours for mean_arr
paths_mean = [Contour.coordinates(curve) for curve in Contour.lines(Contour.levels(contour_result_mean)[1])]
for (x_coords, y_coords) in paths_mean
    # Plot the contour line
    plot!(x_coords, y_coords, color=:black, lw=2, label="")

    # Fill the negative area below the contour line
    plot!(x_coords, y_coords, fillrange=minimum(βh_arr), color=eng_clr, alpha=1.0, lw=0, label="")
end
# hline!(MeanPlot, [0.5], linestyle=:dashdot, color=:red, lw=2.5)
annotate!((-.03, -.03), text("0", :black, 16, :center))

# The white dotted lines
# draw_lines!(MeanPlot, (0,0.4), (0,2), (:x, 0.2), 50; color=:white, markersize=2.)
# draw_lines!(MeanPlot, (0,0.4), (0,2), (:y, 0.4), 50; color=:white, markersize=2.)

println("Saving the plot...")
savefig(plot(MeanPlot, layout=(1, 1), size=(480, 400), left_margin=3mm, bottom_margin=2mm), "RedfieldLoadMeanRate_ω_$(ωe)_δ_$δ.svg") |> println;

# create the plots of two cuts at g = 0.2, and βh = 0.2
g_idx = findfirst(x -> x >= 0.25, g_arr)
βh_idx = findfirst(x -> x >= 0.2, βh_arr)
Pl_g_cut = Pl_arr[g_idx, :]
Pl_βh_cut = Pl_arr[:, βh_idx]

Qh_g_cut = Qh_arr[g_idx, :]
Qh_βh_cut = Qh_arr[:, βh_idx]

Qc_g_cut = Qc_arr[g_idx, :]
Qc_βh_cut = Qc_arr[:, βh_idx]

# plot the g_cut plots and the βh_cut plots separately
g_cut_plot = plot(
    xlim=(βh_arr[begin], βh_arr[end]),
    ylim=(
        -maximum(
            abs.([Pl_g_cut; Qh_g_cut; Qc_g_cut])
            ), 
        maximum(
            abs.([Pl_g_cut; Qh_g_cut; Qc_g_cut])
            )
    ),
    yformatter = (x -> x * 1e4),
    yticks=range(-.0028, .0028, length=5),
    xlabel=L"\beta_h/\beta_c", 
    size=(500, 400),
    margins=3mm,
    xtickfontsize=16,
    ytickfontsize=16,
    guidefontsize=18,
    legendfontsize=16,
)
# find the first indices where Pl, Qh, Qc are negative
acc_idx = findfirst(x -> x < 0, Pl_g_cut) - 1
heat_idx = findfirst(x -> x < 0, Qh_g_cut) - 1
ref_idx = findfirst(x -> x > 0, Qc_g_cut) - 1

# shade the areas between the indices according to their regime
vspan!([βh_arr[begin], βh_arr[acc_idx]], color=eng_clr, alpha=0.5, label="")
vspan!([βh_arr[acc_idx], βh_arr[heat_idx]], color=acc_clr, alpha=0.5, label="")
vspan!([βh_arr[heat_idx], βh_arr[ref_idx]], color=heat_clr, alpha=0.5, label="")
vspan!([βh_arr[ref_idx], βh_arr[end]], color=ref_clr, alpha=0.5, label="")

plot!(
    g_cut_plot, 
    βh_arr, 
    Pl_g_cut, 
    label=L"\langle P_l \rangle",
    lw=4, 
    color="#148109",
)
plot!(
    g_cut_plot, 
    βh_arr, 
    Qh_g_cut, 
    label=L"Q_h", 
    lw=4,
    color="#d62020",
)
plot!(
    g_cut_plot, 
    βh_arr, 
    Qc_g_cut, 
    label=L"Q_c", 
    lw=4,
    color="#1e6cd0",
)
hline!(g_cut_plot, [0.0], linestyle=:dashdot, color=:black, lw=4, label="")

# Save the g_cut plot
# savefig(g_cut_plot, "RedfieldParameterPlotCut_g_$(g_arr[g_idx])_ωe_$(ωe)_δ_$δ.svg") |> println;


βh_cut_plot = plot(
    xlabel=L"g/\hbar\omega_c", 
    xlim=(g_arr[begin], g_arr[end]),ylim=(
        -maximum(
            abs.([Pl_βh_cut; Qh_βh_cut; Qc_βh_cut])
            ), 
        maximum(
            abs.([Pl_βh_cut; Qh_βh_cut; Qc_βh_cut])
            )
    ), 
    yformatter = (x -> x * 1e4),
    yticks=range(-.0016, .0016, length=5),
    size=(500, 400),
    margins=3mm,
    xtickfontsize=16,
    ytickfontsize=16,
    guidefontsize=18,
    legendfontsize=16,
)

# find the first indices where Pl, Qh are negative
acc_idx = findfirst(x -> x < 0, Pl_βh_cut)
heat_idx = findfirst(x -> x < 0, Qh_βh_cut)

# shade the areas between the indices according to their regime
vspan!([g_arr[begin], g_arr[acc_idx]], color=eng_clr, alpha=0.5, label="")
vspan!([g_arr[acc_idx], g_arr[heat_idx]], color=acc_clr, alpha=0.5, label="")
vspan!([g_arr[heat_idx], g_arr[end]], color=heat_clr, alpha=0.5, label="")

plot!(
    βh_cut_plot, 
    g_arr, 
    Pl_βh_cut, 
    label=L"\langle P_l \rangle", 
    lw=4,
    color="#148109",
)
plot!(
    βh_cut_plot, 
    g_arr, 
    Qh_βh_cut, 
    label=L"Q_h", 
    lw=4,
    color="#d62020",
)
plot!(
    βh_cut_plot, 
    g_arr, 
    Qc_βh_cut, 
    label=L"Q_c", 
    lw=4,
    color="#1e6cd0",
)
hline!(βh_cut_plot, [0.0], linestyle=:dashdot, color=:black, lw=4, label="")

# Save the βh_cut plot
# savefig(βh_cut_plot, "RedfieldParameterPlotCut_βh_$(βh_arr[βh_idx])_ωe_$(ωe)_δ_$δ.svg") |> println;

savefig(plot(g_cut_plot, βh_cut_plot, layout=(2, 1), size=(600, 750), left_margin=7mm, bottom_margin=3mm), "RedfieldParameterPlotCuts_ωe_$(ωe)_δ_$δ.svg") |> println;
