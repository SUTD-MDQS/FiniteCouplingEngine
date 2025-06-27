# plot the number mean rate with offset for a given temperature
using JLD2, Plots, Plots.PlotMeasures, LaTeXStrings, Distributed, IterTools
addprocs(40)
@everywhere include("/home/gautham/ThreeLevelLadderEngine/RedfieldParameterSweep/RedfieldThreeLevelLadder.jl")

# the sweep parameters
β_arr = [0.4, 0.45, 0.5, 0.55, 0.6]

# the offset values
offset_arr = 0:20:100

# load the mean rates for all the temperatures
mean_arr = []
scatter_regime_arr = []

ω = 3.0
g = 0.05
δ = 0.0

# Run a simulation if the file is missing
@time @sync @distributed for (βh, offset) in collect(IterTools.product(β_arr, offset_arr))
    filename = "/home/gautham/ThreeLevelLadderEngine/RedfieldOffset/TimeEvolData/OffsetEngine_ω_$(ω)_β_$(βh)_g_$(g)_δ_$(δ)_off_$(offset).jld2"
    if !isfile(filename)
        println("Running simulation for βh = $βh, offset = $offset")
        ergotropy_rate, energy_rate, means, variances, Qh, Qc, t = main(
            βh, ω, g, δ, offset; 
            returnSol=false, 
            steadyStateDynamics=true, 
            savePlots=true, 
            savePlotsAt="/home/gautham/ThreeLevelLadderEngine/RedfieldOffset/Figures/"
        )
        # save the data
        @save filename t means Qh Qc
    end
end

# Load all the simulation data
for βh in β_arr
    mean_offset = []
    scatter_regime = []

    for offset in offset_arr
        filename = "/home/gautham/ThreeLevelLadderEngine/RedfieldOffset/TimeEvolData/OffsetEngine_ω_$(ω)_β_$(βh)_g_$(g)_δ_$(δ)_off_$(offset).jld2"
        println("Loading data for βh = $βh, offset = $offset")
        local t, means, Qh, Qc
        @load filename t means Qh Qc
        push!(mean_offset, linreg(t[end-10:end], means[end-10:end])[2])

        if Qh > 0 && last(mean_offset) > 0 && Qc < 0
            push!(scatter_regime, "engine")
        elseif Qh > 0 && last(mean_offset) < 0 && Qc < 0
            push!(scatter_regime, "accelerator")
        elseif Qh < 0 && last(mean_offset) < 0 && Qc < 0
            push!(scatter_regime, "heater")
        elseif Qh < 0 && last(mean_offset) < 0 && Qc > 0
            push!(scatter_regime, "refrigerator")
        else
            push!(scatter_regime, "unknown")
        end
    end

    push!(mean_arr, mean_offset)
    push!(scatter_regime_arr, scatter_regime)
end


# plot the number mean rate with offset for a given temperature
β_label = reshape([L"β_h:"*"$β" for β in β_arr], (1,:))
# The color palette of our scatter plot
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

p = plot(
    xlabel = L"n_0", 
    ylabel = L"v", 
    xticks = offset_arr,
    yformatter=(val -> round(val * 1e4, digits=1)), 
    ylims=(-0.00034, 0.0002),
    legend = :outerbottom, legendcolumns=3, 
    top_margin=7mm, 
    left_margin=2mm, 
    size=(450, 350),
    xtickfontsize=10,
    ytickfontsize=10,
    guidefontsize=12,
    framestyle = :box,
);
plot!(
    p, 
    offset_arr, 
    mean_arr, 
    label=β_label, 
);
hline!([0], linestyle=:dash, color=:black, label=L"v = 0", linewidth=1.5)
for (mean_offset, scatter_regime) in zip(mean_arr, scatter_regime_arr)
    scatter_colors = [regime_colors[regime] for regime in scatter_regime]  # Map regimes to colors
    scatter_markers = [regime_markers[regime] for regime in scatter_regime]  # Map regimes to markers  

    scatter!(
        p,
        offset_arr, mean_offset,
        color=scatter_colors,
        marker=scatter_markers,
        label="",
        markersize=5
    )
end
annotate!(
    p, 
    3, 0.000255, 
    text(L"\times 10^{-4}", :center, 11, :black)
)

savefig(p, "OffsetPlot.svg")
