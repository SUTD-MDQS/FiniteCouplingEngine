using Plots, JLD, LaTeXStrings, Dierckx, Plots.PlotMeasures
println("Packages loaded successfully")

# Initialize sets to store unique data
ωl_set = Set{Float64}()
β_set = Set{Float64}()

# Read the files in the directory
files = readdir("TimeEvolData/")

# Function to extract β and ω values from the filename
function extract_parameters(filename::String)
    # the data is stored in the following format sol_ω_<>_β_<>.jld
    regex = r"sol_ω_(\d+\.\d+)_β_(\d+\.\d+)\.jld"
    if occursin(regex, filename)
        captures = match(regex, filename).captures
        ω = parse(Float64, captures[1])
        β = parse(Float64, captures[2])
        return (ω, β)
    else
        return nothing
    end
end

# Iterate over the files and load the metadata
for file in files
    # Extract β and ω values from the filename
    params = extract_parameters(file)
    if params !== nothing
        push!(ωl_set, params[1])
        push!(β_set, params[2])
    end
end
println("Files imported successfully!")

# Convert sets to arrays for further processing
ωl_arr = collect(ωl_set) |> sort
β_arr = collect(β_set) |> sort

# the matrices to store the ergotropy rate, ergotropy, energy rate, and whether the engine has population inversion
ergo_rate_mat = zeros(length(ωl_arr), length(β_arr))
carnot_relative_efficiency_mat = zeros(length(ωl_arr), length(β_arr))
energy_rate_mat = zeros(length(ωl_arr), length(β_arr))
refrigerator_efficiency_mat = zeros(length(ωl_arr), length(β_arr))
pop_inversion_mat = zeros(length(ωl_arr), length(β_arr))

println("Computing the matrices from the data...")
# load the data for each value of ωl and β and fill in the matrix values one data at a time
for (i, j) in Iterators.product(eachindex(ωl_arr), eachindex(β_arr))
    ω = ωl_arr[i]
    β = β_arr[j]
	# load the data if it exists
	if !in(("sol_ω_$(ω)_β_$(β).jld"), files)
		println("Data for ω = $ω and β = $β does not exist")
		continue
	end
	
	data = JLD.load("TimeEvolData/sol_ω_$(ω)_β_$(β).jld")
	ergotropy_rate = data["ergotropy_rate"]
	energy_rate = data["energy_rate"]
	pop_inversion = data["pop_inversion"]
    carnot_relative_efficiency = data["carnot_relative_efficiency"]
    refrigerator_efficiency = data["refrigerator_efficiency"]

	ergo_rate_mat[i, j] = ergotropy_rate
	carnot_relative_efficiency_mat[i, j] = carnot_relative_efficiency
	energy_rate_mat[i, j] = energy_rate
	refrigerator_efficiency_mat[i, j] = refrigerator_efficiency
	pop_inversion_mat[i, j] = pop_inversion
end
println("Data computed successfully!")

# # interpolate the energy_rate 2d matrix 
# itp = Spline2D(ωl_arr, β_arr, energy_rate_mat, s=0.01)
x = 0:0.01:maximum(ωl_arr)
y = 0:0.01:maximum(β_arr)
# energy_rate_mat = [itp(i, j) for j in y, i in x]

# # Determine the minimum and maximum values in your ergo_rate_mat and energy_rate_mat
# min_val_ergo_rate, max_val_ergo_rate = minimum(ergo_rate_mat), maximum(ergo_rate_mat)
min_val_energy_rate, max_val_energy_rate = minimum(energy_rate_mat), maximum(energy_rate_mat)

# Calculate the position of zero in the normalized scale
# zero_pos_ergo_rate = (0 - min_val_ergo_rate) / (max_val_ergo_rate - min_val_ergo_rate)
zero_pos_energy_rate = (0 - min_val_energy_rate) / (max_val_energy_rate - min_val_energy_rate)

# Create custom color gradients with white at the calculated zero positions
# custom_rdbu_ergo_rate = cgrad([:blue, :white, :red], [0, zero_pos_ergo_rate, 1])
custom_rdbu_energy_rate = cgrad([:blue, :white, :red], [0, zero_pos_energy_rate, 1])

# # Plot the heatmaps of ergotropy rate, ergotropy, energy rate, and whether the engine has population inversion with respect to ω and β
# p1 = heatmap(ωl_arr, β_arr, ergo_rate_mat, xlabel="ω", ylabel="β", title="Ergotropy Rate", color=custom_rdbu_ergo_rate, legend=:topright, size=(700,600))
# annotate!(p1, 10.5, 2.06, text(L"\dot{W}_e(\rho_{ss})", :center, 12, "black"))

# # plot the (1+ω)β=2 line
# plot!(p1, ωl_arr, 2 ./ (1 .+ ωl_arr), label=L"β_hω_h=β_cω_c", linestyle=:dash, color=:black, xlim=(ωl_arr[1], ωl_arr[end]), ylim=(β_arr[1], β_arr[end]))
# # Plot the curve where pop_inversion_mat changes from true to false
# plot!(p1, ωl_arr, β_false_values, label="Population inversion", linestyle=:solid, color=:green)

# savefig(p1, "Ergotropy rate phase space.svg")
p2 = heatmap(
    ωl_arr, β_arr, energy_rate_mat', 
    xlabel=L"ω_e", ylabel=L"β_h", 
    title="", legend=:topright, 
    color=custom_rdbu_energy_rate, size=(450,400), 
    top_margin=5mm,
    right_margin=3mm, 
 )
# plot the (1+ω)β=2 line
plot!(p2, ωl_arr, 2 ./ (1 .+ ωl_arr), label=L"β_hω_h=β_cω_c", linestyle=:dash , color=:black, xlim=(ωl_arr[1], ωl_arr[end]), ylim=(β_arr[1], β_arr[end]))
annotate!(p2, 10.5, 2.08, text(L"P_l", :center, 12, "black"))
savefig(p2, "EnergyRatePhaseSpace.svg")

# p3 = heatmap(ωl_arr, β_arr, carnot_relative_efficiency_mat, xlabel="ω", ylabel="β", title="carnot relative efficiency", color=:plasma)
# p4 = heatmap(ωl_arr, β_arr, refrigerator_efficiency_mat, xlabel="ω", ylabel="β", title="Refrigeration efficiency", color=:plasma)

# p5 = heatmap(ωl_arr, β_arr, pop_inversion_mat, xlabel="ω", ylabel="β", title="Population Inversion", color=:grays)
# plot!(p5, ωl_arr, 2 ./ (1 .+ ωl_arr), label=L"β_hω_h=β_cω_c", linestyle=:dash, color=:orange)
# savefig(plot(p1, p2, p3, p4, p5, layout=(3, 2), size=(800, 1000)), "Heatmaps.svg")
# println("Plot saved!")