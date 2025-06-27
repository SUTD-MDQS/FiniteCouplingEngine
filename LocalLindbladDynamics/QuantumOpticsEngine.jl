# import the necessary packages
using QuantumOptics, Plots, LaTeXStrings, LinearAlgebra, SparseArrays, Printf

# Function to calculate populations for NLevelBasis
function calculate_populations(sol, engine_basis::NLevelBasis, load_basis::FockBasis)
    # Initialize 2D matrices for populations
    engine_populations = zeros(Float64, engine_basis.N, length(sol))
    load_populations = zeros(Float64, load_basis.N + 1, length(sol))

    # Iterate over sol with index
    for (idx, ρ) in enumerate(sol)
        # Engine populations
        for n in 1:engine_basis.N
            ρproj = projector(nlevelstate(engine_basis, n)) ⊗ identityoperator(load_basis)
            engine_populations[n, idx] = real(expect(ρproj, ρ))
        end

        # Load populations
        for n in 1:load_basis.N+1
            ρproj = identityoperator(engine_basis) ⊗ projector(fockstate(load_basis, n - 1))
            load_populations[n, idx] = real(expect(ρproj, ρ))
        end
    end

    return engine_populations, load_populations
end

# Function to calculate populations for FockBasis
function calculate_populations(sol, engine_basis::FockBasis, load_basis::FockBasis)
    # Initialize 2D matrices for populations
    engine_populations = zeros(Float64, engine_basis.N + 1, length(sol))
    load_populations = zeros(Float64, load_basis.N + 1, length(sol))

    # Iterate over sol with index
    for (idx, ρ) in enumerate(sol)
        # Engine populations
        for n in 1:engine_basis.N+1
            ρproj = projector(fockstate(engine_basis, n - 1)) ⊗ identityoperator(load_basis)
            engine_populations[n, idx] = real(expect(ρproj, ρ))
        end

        # Load populations
        for n in 1:load_basis.N+1
            ρproj = identityoperator(engine_basis) ⊗ projector(fockstate(load_basis, n - 1))
            load_populations[n, idx] = real(expect(ρproj, ρ))
        end
    end

    return engine_populations, load_populations
end


# Compute the entropy of the reduced density matrix
function entropy(ρ)
    # get the eigenvalues of the density matrix
    eigenvalues = eigvals(ρ)

    # filter out the zero eigenvalues to avoid logarithmic divergences
    eigenvalues = eigenvalues[real.(eigenvalues) .> 0]

    # return the entropy
    return real(-sum(eigenvalues .* log.(eigenvalues)))
end

# function to compute the ergotropy of the load
function ergotropy(ρ, H)
    ρ, H = ρ.data, H.data

    # convert ρ and H to dense matrices if they are sparse
    ρ, H = Matrix(ρ), Matrix(H)

    # the energy of the system in its current state
    Ecurrent = real(tr(ρ * H))

    # the density matrix with the eigenvalues sorted in decreasing order
    λ_sorted = sort(eigvals(ρ), by=abs, rev=true)
    ε_sorted = sort(eigvals(H), by=abs)

    # the maximum energy that can be extracted from the system
    Epassive = sum(ε_sorted .* λ_sorted)

    # the ergotropy is the difference between the maximum energy that can be extracted from the system
    E_extractable = Ecurrent - Epassive

    return E_extractable |> real
end

# function to find the heat currents from the hot or cold bath
function heat_currents(sol, H, Jh, Jc, rates)
    # Initialize the heat currents
    Qh = zeros(Float64, length(sol))
    Qc = zeros(Float64, length(sol))

    function D(ρ, J, γ)
        return γ * (J * ρ * dagger(J) - 0.5 * (dagger(J) * J * ρ + ρ * dagger(J) * J))
    end

    # Iterate over sol with index
    for (idx, ρ) in enumerate(sol)
        # Compute the heat currents
        Qh[idx] = real(expect(H, D(ρ, Jh, rates[1]) + D(ρ, dagger(Jh), rates[2])))
        Qc[idx] = real(expect(H, D(ρ, Jc, rates[3]) + D(ρ, dagger(Jc), rates[4])))
    end

    return Qh, Qc
end

# Function to calculate the heat currents coming from the Redfield tensors
function heat_currents_redfield(ρ, H, Rh, Rc)
    # Transformation matrix that brings our density matrix to the eigenbasis of H
    evals, transf_mat = eigen(dense(H).data)
    transf_op = DenseOperator(ρ.basis_l, transf_mat)
    inv_transf_op = DenseOperator(ρ.basis_l, inv(transf_mat))

    # Go to the eigenbasis of the density matrix and vectorize it
    ρ_eb = Ket(ρ.basis_l^2, (inv_transf_op * ρ * transf_op).data[:])
    H_eb = Ket(H.basis_l^2, (Diagonal(evals)[:]))
    
    # Calculate the heat currents in the eigenbasis
    Qh = H_eb.data' * Rh.data * ρ_eb.data
    Qc = H_eb.data' * Rc.data * ρ_eb.data

    return real(Qh), real(Qc)
end

# function to get the linear fit of the ergotropy rate function
function linreg(tspan::Vector{Float64}, arr::Vector{Float64})
    # Ensure that tspan and arr have the same length
    if length(tspan) != length(arr)
        throw(DimensionMismatch("tspan and arr must have the same length"))
    end

    # Example linear regression implementation
    A = [ones(length(tspan)) tspan]
    B = arr

    # Ensure that A and B have compatible dimensions
    if size(A, 1) != length(B)
        throw(DimensionMismatch("A and B must have compatible dimensions"))
    end

    # Perform the linear regression
    coeffs = A \ B
    return coeffs
end

# function to return the mean of an array
function meanN(p::Vector{T}) where T<:Real
    return sum([i * p[i] for i in eachindex(p)])
end

# function to return the variance of an array
function varN(p::Vector{T}) where T<:Real
    m = meanN(p)
    return sum([i^2 * p[i] for i in eachindex(p)]) - m^2
end

# function to find the nonzero elements of a sparse operator with its corresponding basis states
function display_nonzero_elements(operator, threshold)
    function format_complex(x::Complex)
        real_part = real(x)
        imag_part = imag(x)
        formatted_real = @sprintf("%.2e", real_part)
        formatted_imag = @sprintf("%.2e", abs(imag_part))
        sign_imag = imag_part >= 0 ? "+" : "-"
        return "$formatted_real $sign_imag i $formatted_imag"
    end

    sparse_matrix = sparse(operator.data)
    basis_l = operator.basis_l
    basis_r = operator.basis_r
    dim1 = basis_l.bases[1].N + 1
    dim2 = basis_l.bases[2].N + 1

    # Collect non-zero elements with their indices
    elements = [(i, j, value) for (i, j, value) in zip(findnz(sparse_matrix)...) if abs(value) > threshold]

    # Separate diagonal and off-diagonal elements
    diagonal_elements = [(i, j, value) for (i, j, value) in elements if i == j]
    off_diagonal_elements = [(i, j, value) for (i, j, value) in elements if i != j]

    # Sort elements by absolute value
    sorted_diagonal_elements = sort(diagonal_elements, by = x -> -abs(x[3]))
    sorted_off_diagonal_elements = sort(off_diagonal_elements, by = x -> -abs(x[3]))

    # Print diagonal elements
    for (i, j, value) in sorted_diagonal_elements
        bra_state1 = mod(i-1, dim1)
        bra_state2 = mod(div(i-1, dim1), dim2)
        ket_state1 = mod(j-1, dim1)
        ket_state2 = mod(div(j-1, dim1), dim2)
        println("Diagonal value: $value, |$bra_state1⟩⟨$ket_state1| ⊗ |$bra_state2⟩⟨$ket_state2|")
    end

    # Print off-diagonal elements
    for (i, j, value) in sorted_off_diagonal_elements
        bra_state1 = mod(i-1, dim1)
        bra_state2 = mod(div(i-1, dim1), dim2)
        ket_state1 = mod(j-1, dim1)
        ket_state2 = mod(div(j-1, dim1), dim2)
        println("Off-diagonal value: $value, |$bra_state1⟩⟨$ket_state1| ⊗ |$bra_state2⟩⟨$ket_state2|")
    end
end

println("QuantumOpticsEngine.jl ran successfully")

#### TESTING AND DEBUGGING ####
# el, l = [2, 40]

# # define the fock basis
# engineBasis = FockBasis(el)
# loadBasis = FockBasis(l)

# ρ0 = dm(fockstate(engineBasis, 0) ⊗ fockstate(loadBasis, 0))

# display_nonzero_elements(ρ0, 1e-10)