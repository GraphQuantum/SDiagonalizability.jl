# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    minimize_s_bandwidth(g::AbstractGraph, S) -> SBandMinimizationResult
    minimize_s_bandwidth(L::AbstractMatrix{<:Integer}, S) -> SBandMinimizationResult

[TODO: Write here]
"""
function minimize_s_bandwidth(g::AbstractGraph, S::Tuple{Vararg{Integer}})
    _assert_graph_has_defined_s_bandwidth(g)
    res = minimize_s_bandwidth(laplacian_matrix(g), S)
    return SBandMinimizationResult(copy(g), res.S, res.diagonalization, res.band)
end

function minimize_s_bandwidth(L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Integer}})
    _assert_matrix_is_undirected_laplacian(L)

    spec = laplacian_s_spectra(L, S)

    if !spec.s_diagonalizable
        return SBandMinimizationResult(copy(L), S, nothing, Inf)
    end

    pre_eigvecs = hcat(values(spec.s_eigenbases)...)
    # `MatrixBandwidth.jl` uses zero-based indexing, not one-based, so we add 1 here
    band_max = bandwidth(pre_eigvecs'pre_eigvecs) + 1

    k_diag_found = false
    band = 1

    while (!k_diag_found && band < band_max)
        res_temp = _s_spectra_has_bandwidth_at_most_k(spec, band)

        if res_temp.has_band_k_diag
            k_diag_found = true
            diagonalization = res_temp.diagonalization
        else
            band += 1
        end
    end

    if !k_diag_found
        multiplicities = spec.multiplicities
        eigvals = vcat(fill.(keys(multiplicities), values(multiplicities))...)
        diagonalization = Eigen(eigvals, pre_eigvecs)
    end

    return SBandMinimizationResult(copy(L), S, diagonalization, band)
end

"""
    has_s_bandwidth_at_most_k(g::AbstractGraph, S, k) -> SBandRecognitionResult
    has_s_bandwidth_at_most_k(L::AbstractMatrix{<:Integer}, S, k) -> SBandRecognitionResult

[TODO: Write here]
"""
function has_s_bandwidth_at_most_k(g::AbstractGraph, S::Tuple{Vararg{Integer}}, k::Integer)
    _assert_graph_has_defined_s_bandwidth(g)
    res = has_s_bandwidth_at_most_k(laplacian_matrix(g), S, k)
    return SBandRecognitionResult(
        copy(g), res.S, res.diagonalization, k, res.has_band_k_diag
    )
end

function has_s_bandwidth_at_most_k(
    L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Integer}}, k::Integer
)
    _assert_matrix_is_undirected_laplacian(L)
    return _s_spectra_has_bandwidth_at_most_k(laplacian_s_spectra(L, S), k)
end

"""
    is_s_diagonalizable(g::AbstractGraph, S) -> SBandRecognitionResult
    is_s_diagonalizable(L::AbstractMatrix{<:Integer}, S) -> SBandRecognitionResult

[TODO: Write here]
"""
function is_s_diagonalizable(g::AbstractGraph, S::Tuple{Vararg{Integer}})
    _assert_graph_has_defined_s_bandwidth(g)
    res = is_s_diagonalizable(laplacian_matrix(g), S)
    return SBandRecognitionResult(
        copy(g), res.S, res.diagonalization, nv(g), res.has_band_k_diag
    )
end

function is_s_diagonalizable(L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Integer}})
    _assert_matrix_is_undirected_laplacian(L)
    return has_s_bandwidth_at_most_k(L, S, size(L, 1))
end

"""
    _s_spectra_has_bandwidth_at_most_k(spec, k) -> SBandRecognitionResult

[TODO: Write here. Also, comment inline and cite [JP25](@cite)]
"""
function _s_spectra_has_bandwidth_at_most_k(spec::SSpectra, k::Integer)
    L = copy(spec.matrix)
    S = spec.S

    if !spec.s_diagonalizable
        return SBandRecognitionResult(L, S, nothing, k, false)
    end

    n = size(L, 1)

    multiplicities = spec.multiplicities
    s_eigenbases = spec.s_eigenbases

    eigvals = vcat(fill.(keys(multiplicities), values(multiplicities))...)

    # This follows from Johnston and Plosker (2025, p. 320)
    if S == (-1, 1) && n % 4 == 2 && classify_laplacian(L) isa CompleteGraphLaplacian
        if k < n - 1
            res = SBandRecognitionResult(L, S, nothing, k, false)
        else
            eigvecs = hcat(values(s_eigenbases)...)
            diagonalization = Eigen(eigvals, eigvecs)
            res = SBandRecognitionResult(L, S, diagonalization, k, true)
        end

        return res
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
        diagonalization = Eigen(eigvals, eigvecs)
        has_band_k_diag = true
    else
        diagonalization = nothing
        has_band_k_diag = false
    end

    return SBandRecognitionResult(L, S, diagonalization, k, has_band_k_diag)
end
