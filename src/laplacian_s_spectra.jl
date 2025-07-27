# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    laplacian_s_spectra(L, S) -> SSpectra

[TODO: Write here. Also, comment inline]
"""
function laplacian_s_spectra(L::AbstractMatrix{<:Integer}, S::Tuple)
    _assert_matrix_is_undirected_laplacian(L)

    if S == (-1, 0, 1)
        spec = _laplacian_01neg_spectra(L)
    elseif S == (-1, 1)
        spec = _laplacian_1neg_spectra(L)
    else
        throw(ArgumentError("Unsupported entry set S: $S"))
    end

    return spec
end

"""
    check_spectrum_integrality(A) -> SpectrumIntegralResult

Check whether the eigenvalues of `A` are integers (up to floating-point error).

If the eigenvalues are indeed all integers, then an eigenvalue-multiplicity map is
constructed as well.

# Arguments
- `A::AbstractMatrix{<:Integer}`: the matrix whose eigenvalues to check for integrality.

# Returns
- `::SpectrumIntegralResult`: a struct containing the following fields:
    - `matrix`: the `A` matrix, copied to avoid shared mutability.
    - `spectrum_integral`: whether the eigenvalues of `A` are integers.
    - `multiplicities`: a map from each eigenvalue to its multiplicity, sorted first by
        ascending multiplicity then by ascending eigenvalue. (This field is `nothing` if the
        eigenvalues are not all integers.)

# Examples
Confirm that the rotation matrix by ``π/2`` radians counterclockwise is not spectrum
integral (rather, it has eigenvalues ``±i`` [Joy15; p. 1](@cite)):
```jldoctest
julia> R = Int8.([0 -1; 1 0])
2×2 Matrix{Int8}:
 0  -1
 1   0

julia> res = SDiagonalizability.check_spectrum_integrality(R);

julia> res.matrix
2×2 Matrix{Int64}:
 0  -1
 1   0

julia> res.spectrum_integral
false

julia> isnothing(res.multiplicities)
true
```

Confirm that the adjacency matrix of the Petersen graph is spectrum integral, with correct
eigenvalues and multiplicities of ``\\{3: 1, -2: 4, 1: 5\\}`` [Fox09; p. 2](@cite):
```jldoctest
julia> using Graphs

julia> G = smallgraph(:petersen)
{10, 15} undirected simple Int64 graph

julia> A = adjacency_matrix(G)
10×10 SparseArrays.SparseMatrixCSC{Int64, Int64} with 30 stored entries:
 ⋅  1  ⋅  ⋅  1  1  ⋅  ⋅  ⋅  ⋅
 1  ⋅  1  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  1  ⋅  1  ⋅  ⋅  ⋅  1  ⋅  ⋅
 ⋅  ⋅  1  ⋅  1  ⋅  ⋅  ⋅  1  ⋅
 1  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  1
 1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  1  ⋅
 ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  1
 ⋅  ⋅  1  ⋅  ⋅  1  ⋅  ⋅  ⋅  1
 ⋅  ⋅  ⋅  1  ⋅  1  1  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  1  ⋅  1  1  ⋅  ⋅

julia> res = SDiagonalizability.check_spectrum_integrality(A);

julia> res.matrix
10×10 Matrix{Int64}:
 0  1  0  0  1  1  0  0  0  0
 1  0  1  0  0  0  1  0  0  0
 0  1  0  1  0  0  0  1  0  0
 0  0  1  0  1  0  0  0  1  0
 1  0  0  1  0  0  0  0  0  1
 1  0  0  0  0  0  0  1  1  0
 0  1  0  0  0  0  0  0  1  1
 0  0  1  0  0  1  0  0  0  1
 0  0  0  1  0  1  1  0  0  0
 0  0  0  0  1  0  1  1  0  0

julia> res.spectrum_integral
true

julia> res.multiplicities
OrderedCollections.OrderedDict{Int64, Int64} with 3 entries:
  3  => 1
  -2 => 4
  1  => 5
```

# Notes
If an undirected graph with integer edge weights is ``\\{-1, 0, 1\\}``-diagonalizable (or,
more restrictively, ``\\{-1, 1\\}``-diagonalizable), then its Laplacian matrix has integer
eigenvalues [JP25; p. 312](@cite). Hence, validating Laplacian integrality serves as a
useful screening step in this package's principal *S*-bandwidth minimization algorithm.
"""
function check_spectrum_integrality(A::AbstractMatrix{<:Integer})
    A_copy = Matrix{Int}(A) # Avoid shared mutability and cast to `Matrix{Int}`

    eigvals_float = eigvals(A_copy)
    eigvals_int = Int.(round.(real.(eigvals_float)))
    spectrum_integral = isapprox(eigvals_float, eigvals_int)

    # Sort first by ascending multiplicity then by ascending eigenvalue
    if spectrum_integral
        multiplicities = OrderedDict(sort!(collect(counter(eigvals_int)); by=reverse))
    else
        multiplicities = nothing
    end

    return SpectrumIntegralResult(A_copy, spectrum_integral, multiplicities)
end

"""
    _laplacian_01neg_spectra(L) -> SSpectra

[TODO: Write here]
"""
function _laplacian_01neg_spectra(L::AbstractMatrix{<:Integer})
    return _classified_laplacian_01neg_spectra(classify_laplacian(L))
end

"""
    _laplacian_1neg_spectra(L) -> SSpectra
    _laplacian_1neg_spectra(spec) -> SSpectra

[TODO: Write here]
"""
function _laplacian_1neg_spectra(L::AbstractMatrix{<:Integer})
    return _classified_laplacian_1neg_spectra(classify_laplacian(L))
end

function _laplacian_1neg_spectra(spec::SSpectra)
    if spec.S != (-1, 0, 1)
        throw(
            ArgumentError(
                "Expected `{-1, 0, 1}`-spectra` to compute `{-1, 1}`-spectra, got $(spec.S)-spectra",
            ),
        )
    end

    # TODO: Write here. Remember to copy to avoid shared mutability. Fold
    # `_find_indices_1neg` logic from old code into this function.
    return nothing
end

"""
    _classified_laplacian_01neg_spectra(CL) -> SSpectra

[TODO: Write here. Also, comment inline]
"""
function _classified_laplacian_01neg_spectra(CL::ClassifiedLaplacian)
    throw(
        NotImplementedError(
            _classified_laplacian_01neg_spectra, typeof(CL), ClassifiedLaplacian
        ),
    )
end

function _classified_laplacian_01neg_spectra(CL::NullGraphLaplacian)
    return SSpectra(
        CL.matrix,
        (-1, 0, 1),
        OrderedDict{Int,Int}(),
        OrderedDict{Int,Vector{Int}}(),
        OrderedDict{Int,Vector{Int}}(),
        true,
    )
end

function _classified_laplacian_01neg_spectra(CL::EmptyGraphLaplacian)
    L = CL.matrix
    n = size(L, 1)
    S = (-1, 0, 1)

    kernel = stack(pot_kernel_s_eigvecs(n, S))
    # TODO: Justify why this always spans the kernel?
    null_basis = _extract_independent_cols(kernel)

    multiplicities = OrderedDict(0 => n)
    s_eigenspaces = OrderedDict(0 => kernel)
    s_eigenbases = OrderedDict(0 => null_basis)
    s_diagonalizable = true

    return SSpectra(L, S, multiplicities, s_eigenspaces, s_eigenbases, s_diagonalizable)
end

function _classified_laplacian_01neg_spectra(CL::CompleteGraphLaplacian)
    L = CL.matrix
    n = size(L, 1)
    S = (-1, 0, 1)

    #= For every `n ≥ 2`, `Kₙ` (assuming a uniform edge weight of `w`) has eigenvalues 0
    with multiplicity 1 and `n * w` with multiplicity `n - 1`. (`n` is guaranteed to be at
    least 2 in this context due to the Laplacian type casting logic.) =#
    eigenval_nonzero = n * CL.weight

    kernel = ones(Int, n, 1) # The all-ones vector spans the kernel
    null_basis = copy(kernel)

    #= The non-kernel eigenspace necessarily contains all vectors orthogonal to the kernel,
    so we need not filter the output of the generator. =#
    nonkernel_space = stack(pot_nonkernel_s_eigvecs(n, S))
    #= The complete graph is always {-1, 0, 1}-diagonalizable (Johnston and Plosker 2025, p.
    320), so we know that this basis indeed spans the eigenspace. =#
    nonnull_basis = _extract_independent_cols(nonkernel_space)

    multiplicities = OrderedDict(0 => 1, eigenval_nonzero => n - 1)
    s_eigenspaces = OrderedDict(0 => kernel, eigenval_nonzero => nonkernel_space)
    s_eigenbases = OrderedDict(0 => null_basis, eigenval_nonzero => nonnull_basis)
    s_diagonalizable = true

    return SSpectra(L, S, multiplicities, s_eigenspaces, s_eigenbases, s_diagonalizable)
end

function _classified_laplacian_01neg_spectra(CL::ArbitraryGraphLaplacian)
    # TODO: Write here
    return nothing
end

"""
    _classified_laplacian_1neg_spectra(CL) -> SSpectra

[TODO: Write here. Also, comment inline]
"""
function _classified_laplacian_1neg_spectra(CL::ClassifiedLaplacian)
    throw(
        NotImplementedError(
            _classified_laplacian_1neg_spectra, typeof(CL), ClassifiedLaplacian
        ),
    )
end

function _classified_laplacian_1neg_spectra(CL::NullGraphLaplacian)
    return SSpectra(
        CL.matrix,
        (-1, 1),
        OrderedDict{Int,Int}(),
        OrderedDict{Int,Vector{Int}}(),
        OrderedDict{Int,Vector{Int}}(),
        true,
    )
end

function _classified_laplacian_1neg_spectra(CL::EmptyGraphLaplacian)
    L = CL.matrix
    S = (-1, 1)
    n = size(L, 1)

    kernel = stack(pot_kernel_s_eigvecs(n, S))

    #= There exists a {-1, 1}-basis of ℝⁿ if and only if `n = 1` or `n` is even (Johnston
    and Plosker 2025, p. 319). =#
    if n == 1 || iseven(n)
        null_basis = _extract_independent_cols(kernel)
        s_diagonalizable = true
    else
        null_basis = nothing
        s_diagonalizable = false
    end

    multiplicities = OrderedDict(0 => n)
    s_eigenspaces = OrderedDict(0 => kernel)
    s_eigenbases = OrderedDict(0 => null_basis)

    return SSpectra(L, S, multiplicities, s_eigenspaces, s_eigenbases, s_diagonalizable)
end

function _classified_laplacian_1neg_spectra(CL::CompleteGraphLaplacian)
    L = CL.matrix
    n = size(L, 1)
    S = (-1, 1)

    #= For every `n ≥ 2`, `Kₙ` (assuming a uniform edge weight of `w`) has eigenvalues 0
    with multiplicity 1 and `n * w` with multiplicity `n - 1`. (`n` is guaranteed to be at
    least 2 in this context due to the Laplacian type casting logic.) =#
    eigenval_nonzero = n * CL.weight

    kernel = ones(Int, n, 1) # The all-ones vector spans the kernel
    null_basis = copy(kernel)

    #= The non-kernel eigenspace necessarily contains all vectors orthogonal to the kernel,
    so we need not filter the output of the generator. =#
    nonkernel_space = stack(pot_nonkernel_s_eigvecs(n, S))

    #= The complete graph is {-1, 1}-diagonalizable if and only if `n = 1` or `n` is even
    (Johnston and Plosker 2025, p. 320). =#
    if n == 1 || iseven(n)
        nonnull_basis = _extract_independent_cols(nonkernel_space)
        s_diagonalizable = true
    else
        nonnull_basis = nothing
        s_diagonalizable = false
    end

    multiplicities = OrderedDict(0 => 1, eigenval_nonzero => n - 1)
    s_eigenspaces = OrderedDict(0 => kernel, eigenval_nonzero => nonkernel_space)
    s_eigenbases = OrderedDict(0 => null_basis, eigenval_nonzero => nonnull_basis)

    return SSpectra(L, S, multiplicities, s_eigenspaces, s_eigenbases, s_diagonalizable)
end

function _classified_laplacian_1neg_spectra(CL::ArbitraryGraphLaplacian)
    # TODO: Write here
    return nothing
end

"""
    _extract_independent_cols(A) -> Matrix{Int}

Return a (not necessarily unique) independent spanning subset of the columns of `A`.

Computing a rank-revealing (pivoted) QR decomposition of `A`, the scaling coefficients from
the orthogonalization process are used to determine the rank (rather than recompute it with
an SVD), while the pivots are used to extract a spanning set of independent columns.

The rank-revealing Businger–Golub QR algorithm is used for the pivoting strategy, appending
the "most independent" column with respect to the current set of pivots at each step via
Householder transformations [BG65; pp. 269--70](@cite).

# Arguments
- `A::AbstractMatrix{T<:Integer}`: the matrix whose independent columns to extract.

# Returns
- `::Matrix{Int}`: a spanning set of independent columns of `A`.

# Examples
Observe how columns with greater Euclidean norms are given priority in the pivot ordering:
```jldoctest
julia> A = [3  0  0  0  2  1   5   0
            0  3  0  0  2  1  -5   0
            0  0  3  0  2  1   5   4
            0  0  0  3  2  1   0  -4
            0  0  0  0  0  0   0   0]
5×8 Matrix{Int64}:
 3  0  0  0  2  1   5   0
 0  3  0  0  2  1  -5   0
 0  0  3  0  2  1   5   4
 0  0  0  3  2  1   0  -4
 0  0  0  0  0  0   0   0

julia> SDiagonalizability._extract_independent_cols(A)
5×4 Matrix{Int64}:
  5   0  2  3
 -5   0  2  0
  5   4  2  0
  0  -4  2  0
  0   0  0  0
```

# Notes
Since we already need a pivoted QR decomposition to identify independent columns of `A` (or,
rather, to order the columns in such a way that the first `rank(A)` ones are guaranteed to
be independent), it makes sense to use data from the resulting factorization object to
compute the rank of `A` rather than compute a separate SVD. We thus count the nonzero scaling
coefficients—that is, the diagonal entries of the `R` matrix in `A = QR`—to determine the
rank, similarly to how we count the nonzero singular values in an SVD.

It is worth noting that we manually specify a higher relative tolerance for this rank
computation. Further discussion can be found in the [`_rank_rtol`](@ref) documentation, but
in short, a critical part of the formula for `LinearAlgebra.rank`'s default `rtol`
uses the minimum dimension of the input matrix. This may result in rank overestimation for
tall-and-skinny and short-and-fat matrices (precisely the type we expect to encounter when
dealing with all ``\\{-1, 0, 1\\}``-eigenvectors of a Laplacian matrix, which is the
intended use case of this helper function in this package). Our replacement tolerance, on
the other hand, is a widely accepted standard in numerical analysis which uses the maximum
dimension instead [PTVF07; p. 795](@cite).
"""
function _extract_independent_cols(A::AbstractMatrix{<:Integer})
    F = qr(A, ColumnNorm())
    rtol = _rank_rtol(A) # Use a higher tolerance (NumPy's/MATLAB's) than Julia's default

    #= In Julia 1.12+, `LinearAlgebra.rank` dispatches to a method that re-uses an existing
    QR decomposition. For compatibility with v1.10–1.11, we manually define it ourselves in
    `src/utils.jl`. =#
    r = rank(F; rtol=rtol)
    pivots = F.p[1:r] # The first `rank(A)` pivots correspond to independent columns of `A`

    return Matrix{Int}(A[:, pivots]) # Copy to avoid shared mutability
end
