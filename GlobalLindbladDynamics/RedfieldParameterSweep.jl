include("RedfieldThreeLevelEngine.jl")


savePlotsAt = "/home/gautham/ThreeLevelEngineWithQuantumLoad/FiguresSimulation/RedfieldParameterSweep/"
saveTimeEvolDataAt = "/home/gautham/ThreeLevelEngineWithQuantumLoad/TimeEvolData/RedfieldParameterSweep/"


############### Plotting the mean and the variance rates of the load as a function of g ###############
using JLD2, Plots

mean_arr = []
var_arr = []
regimes_arr = []

βh_arr = range(0.01, 1.0, 40)
δ, ωl = 0.0, 3.0
g_arr = range(0.0, 0.75, 40)

# δ_arr = range(-0.2, 0.2, 5)
# δ_arr = [-0.5, -0.25, 0, 0.25, 0.5]

# Iterate over the coupling strength if the mean-variance data does not exist
for βh in βh_arr
    means_g = []
    vars_g = []
    regimes_g = []

    for g in g_arr
        # Set the filename to save the data
        filename = "$(saveTimeEvolDataAt)RedfieldThreeLevelEngine_β_$(βh)_ω_$(ωl)_g_$(g)_δ_$δ.jld2"

        # Check if we already have the data
        if !isfile(filename)
            println("Data does not exist for β=$βh, ω=$ωl, g=$g, δ=$δ. Running the simulation...")
            local ergotropy_rate, energy_rate, means, variances, Qh, Qc, t = main(βh, ωl, g, δ; steadyStateDynamics=true, savePlots=true, savePlotsAt=savePlotsAt)

            @save filename ergotropy_rate energy_rate means variances Qh Qc t
            println("Ran the simulation for β=$βh, ω=$ωl, g=$g, δ=$δ and saved the data.")
        end

        # Load the data from the jld2 file for the given βh and ωl
        println("Loading the data...")
    
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

    # println("Means: ", means)
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

using Dierckx, ImageFiltering

# Marching Squares Algorithm
function marching_squares(matrix::Union{AbstractMatrix{Bool}, AbstractMatrix{Float64}})
    rows, cols = size(matrix)
    contours = []

    # Iterate over each cell in the grid
    for i in 1:(rows - 1)
        for j in 1:(cols - 1)
            # Get the values at the corners of the cell
            tl, tr, bl, br = matrix[i, j], matrix[i, j + 1], matrix[i + 1, j], matrix[i + 1, j + 1]

            # Determine the case based on the corner values
            case = (tl ? 8 : 0) + (tr ? 4 : 0) + (br ? 2 : 0) + (bl ? 1 : 0)

            # Add line segments based on the case
            if case in (1, 14)
                push!(contours, [(j, i + 0.5), (j + 0.5, i + 1)])
            elseif case in (2, 13)
                push!(contours, [(j + 0.5, i + 1), (j + 1, i + 0.5)])
            elseif case in (4, 11)
                push!(contours, [(j + 1, i + 0.5), (j + 0.5, i)])
            elseif case in (8, 7)
                push!(contours, [(j + 0.5, i), (j, i + 0.5)])
            elseif case in (3, 12)
                push!(contours, [(j, i + 0.5), (j + 1, i + 0.5)])
            elseif case in (6, 9)
                push!(contours, [(j + 0.5, i), (j + 0.5, i + 1)])
            elseif case == 5
                push!(contours, [(j, i + 0.5), (j + 0.5, i)])
                push!(contours, [(j + 0.5, i + 1), (j + 1, i + 0.5)])
            elseif case == 10
                push!(contours, [(j + 0.5, i), (j + 1, i + 0.5)])
                push!(contours, [(j, i + 0.5), (j + 0.5, i + 1)])
            end
        end
    end

    return contours
end

# Interpolate the contour points
function interpolate_contours(contours)
    interpolated_contours = []

    for segment in contours
        # Extract x and y coordinates
        xs = [p[1] for p in segment]
        ys = [p[2] for p in segment]

        # Ensure x values are strictly increasing
        if length(xs) > 1
            if all(xs[1:end-1] .<= xs[2:end])  # Check if xs is non-decreasing
                # Interpolate directly if xs is valid
                spline = Spline1D(xs, ys, k=min(3, length(xs) - 1))
                smooth_x = range(minimum(xs), maximum(xs), length=100)
                smooth_y = [spline(x) for x in smooth_x]
                push!(interpolated_contours, (smooth_x, smooth_y))
            else
                # Sort by x if xs is not strictly increasing
                sorted_indices = sortperm(xs)
                xs_sorted = xs[sorted_indices]
                ys_sorted = ys[sorted_indices]
                spline = Spline1D(xs_sorted, ys_sorted, k=min(3, length(xs_sorted) - 1))
                smooth_x = range(minimum(xs_sorted), maximum(xs_sorted), length=100)
                smooth_y = [spline(x) for x in smooth_x]
                push!(interpolated_contours, (smooth_x, smooth_y))
            end
        else
            # If only one point, just add it as-is
            push!(interpolated_contours, (xs, ys))
        end
    end

    return interpolated_contours
end

# Plot the matrix and interpolated contours
function plot_marching_squares!(fig, matrix::Union{AbstractMatrix{Bool}, AbstractMatrix{Float64}})
    # Get the contours using marching squares
    contours = marching_squares(matrix)

    # Interpolate the contours
    interpolated_contours = interpolate_contours(contours)

    # Plot the heatmap
    # fig = heatmap(matrix, color=:viridis, xlabel="X", ylabel="Y", title="Marching Squares Contours")

    # Overlay the interpolated contours with a single label
    first = true
    for (smooth_x, smooth_y) in interpolated_contours
        if first
            plot!(fig, smooth_x, smooth_y, color=:red, linewidth=3, label="Contour 1")
            first = false
        else
            plot!(fig, smooth_x, smooth_y, color=:red, linewidth=3, label="")
        end
    end
end

function plot_raw_contours!(fig, matrix::AbstractMatrix{Bool})
    # Get the contours using marching squares
    contours = marching_squares(matrix)

    # Plot the heatmap
    # fig = heatmap(matrix, color=:viridis, xlabel="X", ylabel="Y", title="Raw Marching Squares Contours")

    # Overlay the raw contours
    for segment in contours
        xs = [p[1] for p in segment]
        ys = [p[2] for p in segment]
        plot!(fig, xs, ys, color=:blue, linewidth=2, label="")
    end

end


# Plot the mean and variance rates
MeanPlot = heatmap(
    g_arr, 
    βh_arr,
    mean_arr, 
    xlabel=L"g", 
    ylabel=L"\beta_h", 
    # L"\dot{⟨n_l⟩}"
    title=L"\dot{⟨n_l⟩}"*" for ωl, δ= $(ωl), $δ", 
    color=rdbu_mean
)
VarPlot = heatmap(
    g_arr, 
    βh_arr,
    var_arr,
    xlabel=L"g", 
    ylabel=L"\beta_h", 
    # L"⟨Δ\dot{n}_l⟩^2"
    title=L"⟨Δ\dot{n}_l⟩^2"*" for ωl, δ= $(ωl), $δ", 
    color=rdbu_var
)
contour!(
    MeanPlot,
    g_arr, 
    βh_arr,
    mean_arr, 
    levels=[0]
)
contour!(MeanPlot, g_arr, βh_arr, engine_matrix, levels=[0.00005], linewidth=2, color=:black, label="Engine")
contour!(MeanPlot, g_arr, βh_arr, accelerator_matrix, levels=[0.00005], linewidth=2, color=:red, label="Accelerator")
contour!(MeanPlot, g_arr, βh_arr, refrigerator_matrix, levels=[0.00005], linewidth=2, color=:blue, label="Refrigerator")
contour!(MeanPlot, g_arr, βh_arr, heater_matrix, levels=[0.00005], linewidth=2, color=:green, label="Heater")
contour!(MeanPlot, g_arr, βh_arr, imfilter(engine_matrix, Kernel.gaussian(2.5)), levels=[0.00005], linewidth=2, color=:black, label="Engine")

# println("Saving the plot...")
# savefig(MeanPlot, "MeanPlot.svg")
# savefig(VarPlot, "VarPlot.svg")
# savefig(plot(MeanPlot, VarPlot, layout=(2, 1), size=(600, 1000), margin=5mm), "RedfieldLoadMeanVarianceRate_ω_$(ωl)_δ_$δ.svg") |> println;

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
