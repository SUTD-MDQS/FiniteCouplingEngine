#!/bin/bash

# this is a bash script that runs the Sample.jl in parallel using the GNU parallel tool. Note that the figure generated arent combined into a single figure, rather saved as separate figures.
start_time=$(date +%s.%3N)

# we have 80 different values of beta and omega
beta=$(seq 0.1 0.02375 2.0)
omega=$(seq 0.0 0.125 10.0)

nthreads=$(nproc)
nthreads=150

echo "Starting parallel execution of julia scripts"

# loop over the values of a and run the simulation in parallel
parallel -j $nthreads JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 julia QuantumOpticsEngine.jl {1} {2} ::: $beta ::: $omega


end_time=$(date +%s.%3N)

total_time=$(echo "$end_time - $start_time" | bc)

echo "Total time taken: $total_time seconds"
echo "Number of parallel computations performed: $size"
