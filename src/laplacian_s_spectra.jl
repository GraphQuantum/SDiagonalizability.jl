# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    laplacian_s_spectra(L, S) -> SSpectra

[TODO: Write here. Also, comment inline]
"""
function laplacian_s_spectra(L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Integer}})
    _assert_matrix_is_undirected_laplacian(L)

    S = _sort_tuple(S)

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
integral (rather, it has eigenvalues ``±i`` [Joy15, p. 1]):
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
eigenvalues and multiplicities of ``\\{3: 1, -2: 4, 1: 5\\}`` [Fox09, p. 2]:
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
eigenvalues [JP25, p. 312]. Hence, validating Laplacian integrality serves as a useful
screening step in this package's principal *S*-bandwidth minimization algorithm.

# References

- [Fox09](@cite): J. Fox. *Lecture 19: The Petersen graph and Moore graphs*. Lecture notes,
    MAT 307: Combinatorics (2009). Accessed: 2025-07-25.
    https://math.mit.edu/~fox/MAT307.html.
- [Joy15](@cite): D. Joyce. *Rotations and complex eigenvalues*. Lecture notes, Math 130:
    Linear Algebra (2015). http://aleph0.clarku.edu/~ma130/complexeigen.pdf.
- [JP25](@cite): N. Johnston and S. Plosker. *Laplacian {−1,0,1}- and {−1,1}-diagonalizable
    graphs*. Linear Algebra and its Applications **704**, 309–39 (2025).
    https://doi.org/10.1016/j.laa.2024.10.016.
"""
function check_spectrum_integrality(A::AbstractMatrix{<:Integer})
    A_copy = Matrix{Int}(A) # Avoid shared mutability and cast to `Matrix{Int}`
    n = size(A_copy, 1)

    eigvals_float = eigvals(A_copy)
    eigvals_int = Vector{Int}(undef, n)
    spectrum_integral = true
    i = 0

    while (spectrum_integral && i < n)
        eigval_float = eigvals_float[i += 1]
        eigval_int = round(Int, real(eigval_float))
        eigvals_int[i] = eigval_int

        if !isapprox(eigval_float, eigval_int; atol=1e-8, rtol=1e-5)
            spectrum_integral = false
        end
    end

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

[TODO: Write here]
"""
function _laplacian_1neg_spectra(L::AbstractMatrix{<:Integer})
    return _classified_laplacian_1neg_spectra(classify_laplacian(L))
end

"""
    _classified_laplacian_01neg_spectra(CL) -> SSpectra

[TODO: Write here. Also, comment inline]
"""
function _classified_laplacian_01neg_spectra(::T) where {T<:ClassifiedLaplacian}
    throw(NotImplementedError(_classified_laplacian_01neg_spectra, T, ClassifiedLaplacian))
end

function _classified_laplacian_01neg_spectra(CL::NullGraphLaplacian)
    return SSpectra(
        CL.matrix,
        (-1, 0, 1),
        OrderedDict{Int,Int}(),
        OrderedDict{Int,Matrix{Int}}(),
        OrderedDict{Int,Matrix{Int}}(),
        true,
    )
end

function _classified_laplacian_01neg_spectra(CL::EmptyGraphLaplacian)
    L = CL.matrix
    n = size(L, 1)
    S = (-1, 0, 1)

    kernel = stack(pot_kernel_s_eigvecs(n, S))
    null_basis = _extract_independent_cols(kernel) # This will always form a basis of ℝⁿ

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
    nonkernel_vecs = pot_nonkernel_s_eigvecs(n, S)

    if isempty(nonkernel_vecs) # In case `S = (-1, 0, 1)` and `n` is odd
        nonkernel_space = Matrix{Int}(undef, n, 0)
    else
        nonkernel_space = stack(pot_nonkernel_s_eigvecs(n, S))
    end

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
    # The logic for arbitrary graphs is identical for both `{-1, 0, 1}` and `{-1, 1}`
    return _classified_laplacian_s_spectra(CL, (-1, 0, 1))
end

"""
    _classified_laplacian_1neg_spectra(CL) -> SSpectra

[TODO: Write here. Also, comment inline]
"""
function _classified_laplacian_1neg_spectra(::T) where {T<:ClassifiedLaplacian}
    throw(NotImplementedError(_classified_laplacian_1neg_spectra, T, ClassifiedLaplacian))
end

function _classified_laplacian_1neg_spectra(CL::NullGraphLaplacian)
    return SSpectra(
        CL.matrix,
        (-1, 1),
        OrderedDict{Int,Int}(),
        OrderedDict{Int,Matrix{Int}}(),
        OrderedDict{Int,Matrix{Int}}(),
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
    nonkernel_vecs = pot_nonkernel_s_eigvecs(n, S)

    if isempty(nonkernel_vecs) # In case `S = (-1, 1)` and `n` is odd
        nonkernel_space = Matrix{Int}(undef, n, 0)
    else
        nonkernel_space = stack(pot_nonkernel_s_eigvecs(n, S))
    end

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
    # The logic for arbitrary graphs is identical for both `{-1, 0, 1}` and `{-1, 1}`
    return _classified_laplacian_s_spectra(CL, (-1, 1))
end

"""
    _classified_laplacian_s_spectra(CL, S) -> SSpectra

[TODO: Write here]
"""
function _classified_laplacian_s_spectra(
    CL::ArbitraryGraphLaplacian, S::Tuple{Vararg{Integer}}
)
    L = CL.matrix
    S = _sort_tuple(S)
    res = check_spectrum_integrality(L)

    if !res.spectrum_integral
        return SSpectra(
            L,
            S,
            OrderedDict{Int,Int}(),
            OrderedDict{Int,Matrix{Int}}(),
            OrderedDict{Int,Matrix{Int}}(),
            false,
        )
    end

    n = size(L, 1)

    #= `multiplicities` is only `nothing` when `spectrum_integral` is false, so we
    explicitly assert its type here for compiler inference during static analysis. =#
    multiplicities = res.multiplicities
    @assert !isnothing(multiplicities)

    eigvals_nonzero = filter(!iszero, keys(multiplicities))
    s_eigenspaces = OrderedDict{Int,AbstractMatrix{Int}}(
        # Initialize resizeable matrices to store an unknown number of eigenvectors
        eigval => ElasticMatrix{Int}(undef, n, 0) for eigval in eigvals_nonzero
    )
    s_eigenbases = OrderedDict{Int,Union{Nothing,Matrix{Int}}}()

    # The kernel is simpler and handled differently from the other eigenspaces
    if multiplicities[0] == 1 # The all-ones vector spans the kernel
        kernel = ones(Int, n, 1)
    else # Check all `S`-vectors in ℝⁿ, unique up to span
        kernel_vecs = Iterators.filter(v -> iszero(L * v), pot_kernel_s_eigvecs(n, S))

        if isempty(kernel_vecs)
            kernel = Matrix{Int}(undef, n, 0)
        else
            kernel = stack(kernel_vecs)
        end
    end

    s_eigenspaces[0] = kernel

    for v in pot_nonkernel_s_eigvecs(n, S)
        for eigval in eigvals_nonzero
            if L * v == eigval * v
                append!(s_eigenspaces[eigval], v)
                break
            end
        end
    end

    for eigval in keys(multiplicities)
        eigenspace = s_eigenspaces[eigval]

        if isempty(eigenspace)
            eigenbasis = nothing
        else
            pot_eigenbasis = _extract_independent_cols(eigenspace)

            if size(pot_eigenbasis, 2) < multiplicities[eigval]
                eigenbasis = nothing
            else
                eigenbasis = pot_eigenbasis
            end
        end

        s_eigenbases[eigval] = eigenbasis
    end

    s_diagonalizable = all(!isnothing, values(s_eigenbases))

    return SSpectra(L, S, multiplicities, s_eigenspaces, s_eigenbases, s_diagonalizable)
end
