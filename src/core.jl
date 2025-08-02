# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    s_bandwidth(g::AbstractGraph, S) -> SBandMinimizationResult
    s_bandwidth(L::AbstractMatrix{<:Integer}, S) -> SBandMinimizationResult

[TODO: Write here. Also, comment inline and cite [JP25](@cite).]
"""
function s_bandwidth(g::AbstractGraph, S::Tuple{Vararg{Integer}})
    _assert_graph_has_defined_s_bandwidth(g)

    res = s_bandwidth(laplacian_matrix(g), S)

    return SBandMinimizationResult(copy(g), res.S, res.s_diagonalization, res.s_bandwidth)
end

function s_bandwidth(L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Integer}})
    _assert_matrix_is_undirected_laplacian(L)

    S = _sort_tuple(S)
    spec = laplacian_s_spectra(L, S)

    if !spec.s_diagonalizable
        return SBandMinimizationResult(copy(L), S, nothing, Inf)
    end

    n = size(L, 1)

    multiplicities = spec.multiplicities
    s_eigenspaces = spec.s_eigenspaces
    s_eigenbases = spec.s_eigenbases

    eigvals = vcat(fill.(keys(multiplicities), values(multiplicities))...)
    eigvecs = Matrix{Int}(undef, n, n)
    idx = 1

    # Using results from Johnston and Plosker (2025, pp. 319–323)
    if classify_laplacian(L) isa CompleteGraphLaplacian
        if S == (-1, 0, 1) && 2 < n <= 20 && n % 4 != 0 # Subcase of Conjecture 2
            k = 2
        elseif S == (-1, 1) && n % 4 == 2 # Theorem 1
            k = n - 1
        else # Back to the general case
            k = 1
        end
    else
        k = 1
    end

    for (eigval, multi) in multiplicities
        k_basis = find_k_orthogonal_basis(
            s_eigenspaces[eigval], multi, k, s_eigenbases[eigval]
        )

        while isnothing(k_basis)
            k += 1
            k_basis = find_k_orthogonal_basis(
                s_eigenspaces[eigval], multi, k, s_eigenbases[eigval]
            )
        end

        eigvecs[:, idx:((idx += multi) - 1)] .= k_basis
    end

    s_diagonalization = Eigen(eigvals, eigvecs)

    return SBandMinimizationResult(copy(L), S, s_diagonalization, k)
end

"""
    has_s_bandwidth_at_most_k(g::AbstractGraph, S, k) -> SBandRecognitionResult
    has_s_bandwidth_at_most_k(L::AbstractMatrix{<:Integer}, S, k) -> SBandRecognitionResult

[TODO: Write here. Also, comment inline and cite [JP25](@cite).]
"""
function has_s_bandwidth_at_most_k(g::AbstractGraph, S::Tuple{Vararg{Integer}}, k::Integer)
    _assert_graph_has_defined_s_bandwidth(g)

    res = has_s_bandwidth_at_most_k(laplacian_matrix(g), S, k)

    return SBandRecognitionResult(
        copy(g), res.S, res.s_diagonalization, k, res.s_band_at_most_k
    )
end

function has_s_bandwidth_at_most_k(
    L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Integer}}, k::Integer
)
    _assert_matrix_is_undirected_laplacian(L)

    S = _sort_tuple(S)
    spec = laplacian_s_spectra(L, S)
    L_copy = copy(spec.matrix)

    if !spec.s_diagonalizable
        return SBandRecognitionResult(L_copy, S, nothing, k, false)
    end

    n = size(L_copy, 1)

    multiplicities = spec.multiplicities
    s_eigenbases = spec.s_eigenbases

    eigvals = vcat(fill.(keys(multiplicities), values(multiplicities))...)

    # Using results from Johnston and Plosker (2025, pp. 319–323)
    if classify_laplacian(L_copy) isa CompleteGraphLaplacian
        if S == (-1, 0, 1) && k == 1 && 2 < n <= 20 && n % 4 != 0 # Subcase of Conjecture 2
            return SBandRecognitionResult(L_copy, S, nothing, k, false)
        elseif S == (-1, 1) && n % 4 == 2 # Theorem 1
            if k < n - 1
                res = SBandRecognitionResult(L_copy, S, nothing, k, false)
            else
                eigvecs = hcat(values(s_eigenbases)...)
                s_diagonalization = Eigen(eigvals, eigvecs)
                res = SBandRecognitionResult(L_copy, S, s_diagonalization, k, true)
            end

            return res
        end
    end

    s_eigenspaces = spec.s_eigenspaces

    eigvecs = Matrix{Int}(undef, n, n)

    res = iterate(multiplicities)
    k_basis_found = true
    idx = 1

    while (k_basis_found && idx <= n)
        ((eigval, multi), state) = res
        k_basis = find_k_orthogonal_basis(
            s_eigenspaces[eigval], multi, k, s_eigenbases[eigval]
        )

        if isnothing(k_basis)
            k_basis_found = false
        else
            eigvecs[:, idx:(idx + multi - 1)] .= k_basis
            res = iterate(multiplicities, state)
            idx += multi
        end
    end

    if k_basis_found
        s_diagonalization = Eigen(eigvals, eigvecs)
        s_band_at_most_k = true
    else
        s_diagonalization = nothing
        s_band_at_most_k = false
    end

    return SBandRecognitionResult(L_copy, S, s_diagonalization, k, s_band_at_most_k)
end

"""
    is_s_diagonalizable(g::AbstractGraph, S) -> SDiagonalizabilityResult
    is_s_diagonalizable(L::AbstractMatrix{<:Integer}, S) -> SDiagonalizabilityResult

[TODO: Write here]
"""
function is_s_diagonalizable(g::AbstractGraph, S::Tuple{Vararg{Integer}})
    _assert_graph_has_defined_s_bandwidth(g)

    res = is_s_diagonalizable(laplacian_matrix(g), S)

    return SDiagonalizabilityResult(
        copy(g), res.S, res.s_diagonalization, res.has_s_diagonalization
    )
end

function is_s_diagonalizable(L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Integer}})
    _assert_matrix_is_undirected_laplacian(L)

    S = _sort_tuple(S)
    spec = laplacian_s_spectra(L, S)
    L_copy = copy(spec.matrix)

    if spec.s_diagonalizable
        multiplicities = spec.multiplicities
        eigvals = vcat(fill.(keys(multiplicities), values(multiplicities))...)
        eigvecs = hcat(values(spec.s_eigenbases)...)
        s_diagonalization = Eigen(eigvals, eigvecs)
        has_s_diagonalization = true
    else
        s_diagonalization = nothing
        has_s_diagonalization = false
    end

    return SDiagonalizabilityResult(L_copy, S, s_diagonalization, has_s_diagonalization)
end
