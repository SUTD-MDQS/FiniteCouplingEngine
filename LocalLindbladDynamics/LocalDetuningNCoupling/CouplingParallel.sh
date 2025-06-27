#!/bin/bash

# This is a bash script that runs the QuantumOpticsEngine.jl in parallel using GNU Parallel.
start_time=$(date +%s.%3N)

# Define the parameters
beta_vals=(0.45)
omega_vals=(1.0)
g_vals=(0.01 0.05 0.1 0.5 1.0)

# Export the parameters to be used by GNU Parallel
export beta_vals
export omega_vals
export g_vals

# Function to run the Julia script with given parameters
run_julia() {
    beta=$1
    omega=$2
    g=$3
    julia -e "include(\"ThreeLevelEngineWithQuantumLoad.jl\"); main($beta, $omega, $g)"
}

export -f run_julia

# Number of threads to use
nthreads=$(nproc)

# Run the simulation in parallel
parallel -j $nthreads OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 run_julia ::: "${beta_vals[@]}" ::: "${omega_vals[@]}" ::: "${g_vals[@]}"

wait

end_time=$(date +%s.%3N)

total_time=$(echo "$end_time - $start_time" | bc)

echo "Total time taken: $total_time seconds"
echo "Number of parallel computations performed: ${#beta_vals[@]} * ${#omega_vals[@]} * ${#g_vals[@]}"
