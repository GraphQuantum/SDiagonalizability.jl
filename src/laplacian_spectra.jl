# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    struct SpectrumIntegralResult

Data on whether a matrix is spectrum integral (i.e., whether its eigenvalues are integers).

This struct also contains a map from each eigenvalue to its multiplicity, provided that the
eigenvalues are indeed all integers. (Otherwise, the associated field is simply `nothing`.)

# Fields
- `matrix::Matrix{Int}`: the matrix whose eigenvalues and their integrality are of interest.
- `spectrum_integral::Bool`: whether the eigenvalues of `matrix` are all integers.
- `multiplicities::Union{Nothing,OrderedDict{Int,Int}}`: a map from each eigenvalue to its
    multiplicity, sorted first by ascending multiplicity then by ascending eigenvalue. (This
    field is `nothing` if `spectrum_integral` is false.)

# Notes
If an undirected graph with integer edge weights is `{-1,0,1}`-diagonalizable (or, more
restrictively, `{-1,1}`-diagonalizable), then its Laplacian matrix has integer eigenvalues
[JP25; p. 312](@cite). Hence, validating Laplacian integrality serves as a useful screening
step in this package's principal *S*-bandwidth minimization algorithm.
"""
struct SpectrumIntegralResult
    matrix::Matrix{Int}
    spectrum_integral::Bool
    multiplicities::Union{Nothing,OrderedDict{Int,Int}}

    function SpectrumIntegralResult(
        matrix::Matrix{Int},
        spectrum_integral::Bool,
        multiplicities::Union{Nothing,OrderedDict{Int,Int}},
    )
        if spectrum_integral && isnothing(multiplicities)
            throw(
                ArgumentError(
                    "`multiplicities` is only `nothing` as a sentinel value " *
                    "when `spectrum_integral` is false.",
                ),
            )
        end

        return new(matrix, spectrum_integral, multiplicities)
    end
end

"""
    struct _Eigenspace01Neg

Data on an eigenspace of some Laplacian matrix with respect to the `{-1,0,1}`-spectrum.

# Fields
- `multiplicity::Int`: the multiplicity of the eigenspace.
- `eigvecs_01neg::AbstractMatrix{Int}`: the `{-1,0,1}`-eigenvectors in this eigenspace,
    stored as columns of a matrix.
- `indices_1neg::Vector{Int}`: the indices of all `{-1,1}`-columns in `eigvecs_01neg`.
- `basis_01neg::Union{Nothing,Matrix{Int}}`: a `{-1,0,1}`-basis for the eigenspace, if one
    exists. (This field is `nothing` if no such basis exists.)
- `basis_1neg::Union{Nothing,Matrix{Int}}`: a `{-1,1}`-basis for the eigenspace, if one
    exists. (This field is `nothing` if no such basis exists.)

# Notes
This helper struct is, above all else, used for the sake of efficient data storage in
user-facing structs, with the [`LaplacianSpectrum01Neg`](@ref) constructor flattening data
stored across several `_Eigenspace01Neg` instances into top-level fields. This justifies the
exclusion of the `eigval` field, which is redundant in the context of
[`LaplacianSpectrum01Neg`](@ref) (the keys of the `eigspaces_01neg::{Int,_Eigenspace01Neg}`
map already specify the eigenvalues).

See the [`_find_indices_1neg`](@ref) documentation if seeking justification for why the
`indices_1neg` field is a `Vector{Int}` rather than a `BitVector`. (In short, very few
columns have entries exclusively from `{-1,1}`, so `Vector{Int}`'s consume less memory.)
"""
struct _Eigenspace01Neg
    multiplicity::Int
    eigvecs_01neg::AbstractMatrix{Int}
    indices_1neg::Vector{Int}
    basis_01neg::Union{Nothing,Matrix{Int}}
    basis_1neg::Union{Nothing,Matrix{Int}}
end

"""
    struct LaplacianSpectrum01Neg{R,S<:R,T<:R,U<:R,V<:R}

Data on the `{-1,0,1}`-spectrum of some Laplacian matrix.

Each instance encapsulates the Laplacian itself, whether the Laplacian is `{-1,0,1}`- and
`{-1,1}`-diagonalizable, each eigenvalue's multiplicity, and the associated `{-1,0,1}`- and
`{-1,1}`-eigenvectors, and eigenbases. If the Laplacian is not spectrum integral, then it
cannot have any `{-1,0,1}`-eigenvectors outside of the kernel [JP25; p. 312](@cite), so much
of this data is simply `nothing`.

# Type Parameters
- `R<:Union{Nothing,OrderedDict{Int}}`: the "common upper bound" for all five map-like
    fields, ensuring that they are either all of the same supertype. Either `Nothing` (if
    and only if the Laplacian is not spectrum integral, since we cannot map non-integer
    eigenvalues to data) or `OrderedDict{Int}` (if all eigenvalues are indeed integers).
- `S<:R`: the type of the `multiplicities` field. Either `Nothing` or an
    `OrderedDict{Int,Int}`.
- `T<:R`: the type of the `eigvecs_01neg` field. Either `Nothing` or an
    `OrderedDict{Int,AbstractMatrix{Int}}`.
- `U<:R`: the type of the `indices_1neg` field. Either `Nothing` or an
    `OrderedDict{Int,Vector{Int}}`.
- `V<:R`: the type of the `bases_01neg` and `bases_1neg` fields. Either `Nothing` or an
    `OrderedDict{Int,Union{Nothing,Matrix{Int}}}`.

# Fields
- `laplacian_matrix::Matrix{Int}`: the Laplacian matrix.
- `diagonalizable_01neg::Bool`: whether the Laplacian is `{-1,0,1}`-diagonalizable.
- `diagonalizable_1neg::Bool`: whether the Laplacian is `{-1,1}`-diagonalizable.
- `multiplicities::S`: a map from each eigenvalue to its multiplicity, sorted first by
    ascending multiplicity then by ascending eigenvalue. (If the Laplacian is not spectrum
    integral, this field is `nothing`.)
- `eigvecs_01neg::T`: a map from each eigenvalue to its `{-1,0,1}`-eigenvectors, stored as
    columns of a matrix. (If the Laplacian is not spectrum integral, this field is
    `nothing`.)
- `indices_1neg::U`: a map from each eigenvalue to the indices of the `{-1,1}`-columns in
    `eigvecs_01neg[eigval]`. (If the Laplacian is not spectrum integral, this field is
    `nothing`.)
- `bases_01neg::V`: a map from each eigenvalue to a `{-1,0,1}`-basis for the corresponding
    eigenspace if one exists, or to `nothing` otherwise. (If the Laplacian is not spectrum
    integral, this field is `nothing`.)
- `bases_1neg::V`: a map from each eigenvalue to a `{-1,1}`-basis for the corresponding
    eigenspace if one exists, or to `nothing` otherwise. (If the Laplacian is not spectrum
    integral, this field is `nothing`.)

# Constructors
- `LaplacianSpectrum01Neg(::Matrix{Int}, ::Bool, ::Bool, ::Union{Nothing,OrderedDict{Int,_Eigenspace01Neg}})`:
    constructs a `LaplacianSpectrum01Neg` instance from the given Laplacian matrix and
    associated data on its `{-1,0,1}`-eigenspaces. When the last argument is not `nothing`
    (indicating spectrum integrality), this constructors flattens the data stored in the
    [`_Eigenspace01Neg`](@ref) instances and exposes it via top-level fields. (See the
    [`_Eigenspace01Neg`](@ref) documentation for more details.)

# Notes
TODO: Discuss potential use cases for this beyond the `s_bandwidth` algorithm
"""
struct LaplacianSpectrum01Neg{R<:Union{Nothing,OrderedDict{Int}},S<:R,T<:R,U<:R,V<:R}
    laplacian_matrix::Matrix{Int}
    diagonalizable_01neg::Bool
    diagonalizable_1neg::Bool
    multiplicities::S
    eigvecs_01neg::T
    indices_1neg::U
    bases_01neg::V
    bases_1neg::V

    function LaplacianSpectrum01Neg(
        laplacian_matrix::Matrix{Int},
        diagonalizable_01neg::Bool,
        diagonalizable_1neg::Bool,
        eigspaces_01neg::Union{Nothing,OrderedDict{Int,_Eigenspace01Neg}},
    )
        #= {-1,1}-diagonalizability is a strict subset of {-1,0,1}-diagonalizability, so any
        matrix that is {-1,1}-diagonalizable must also be {-1,0,1}-diagonalizable. =#
        if diagonalizable_1neg && !diagonalizable_01neg
            throw(
                ArgumentError(
                    "If `diagonalizable_1neg` is `true`, then `diagonalizable_01neg` " *
                    "must be `true` as well",
                ),
            )
        end

        #= Call a helper that dispatches to different methods depending on the type of
        `eigspaces_01neg`. The eigenspace data is flattened therein. =#
        return (_helper_constructor(
            laplacian_matrix, diagonalizable_01neg, diagonalizable_1neg, eigspaces_01neg
        ))
    end

    function _helper_constructor(
        laplacian_matrix::Matrix{Int},
        diagonalizable_01neg::Bool,
        diagonalizable_1neg::Bool,
        ::Nothing,
    )
        #= `eigspaces_01neg` is `nothing` if and only if the Laplacian is not spectrum
        integral, in which case the Laplacian cannot be {-1,0,1}-diagonalizable either. =#
        if diagonalizable_01neg
            throw(
                ArgumentError(
                    "If `eigspaces_01neg` is `nothing`, then `diagonalizable_01neg` must " *
                    " be false",
                ),
            )
        end

        R = S = T = U = V = Nothing

        #= We call `S()`, `T()`, etc. to signify the relationship between the fields and
        their type parameters, but this is really just equivalent to passing `nothing`. =#
        return new{R,S,T,U,V}(
            laplacian_matrix,
            diagonalizable_01neg,
            diagonalizable_1neg,
            S(),
            T(),
            U(),
            V(),
            V(),
        )
    end

    function _helper_constructor(
        laplacian_matrix::Matrix{Int},
        diagonalizable_01neg::Bool,
        diagonalizable_1neg::Bool,
        eigspaces_01neg::OrderedDict{Int,_Eigenspace01Neg},
    )
        R = OrderedDict{Int}
        S = OrderedDict{Int,Int}
        T = OrderedDict{Int,AbstractMatrix{Int}}
        U = OrderedDict{Int,Vector{Int}}
        V = OrderedDict{Int,Union{Nothing,Matrix{Int}}}

        #= Here we flatten the data stored in the `_Eigenspace01Neg` instances to expose it
        via top-level fields of `LaplacianSpectrum01Neg`. =#
        multiplicities = S(
            eigval => eigspace.multiplicity for (eigval, eigspace) in eigspaces_01neg
        )
        eigvecs_01neg = T(
            eigval => eigspace.eigvecs_01neg for (eigval, eigspace) in eigspaces_01neg
        )
        indices_1neg = U(
            eigval => eigspace.indices_1neg for (eigval, eigspace) in eigspaces_01neg
        )
        bases_01neg = V(
            eigval => eigspace.basis_01neg for (eigval, eigspace) in eigspaces_01neg
        )
        bases_1neg = V(
            eigval => eigspace.basis_1neg for (eigval, eigspace) in eigspaces_01neg
        )

        return new{R,S,T,U,V}(
            laplacian_matrix,
            diagonalizable_01neg,
            diagonalizable_1neg,
            multiplicities,
            eigvecs_01neg,
            indices_1neg,
            bases_01neg,
            bases_1neg,
        )
    end
end

"""
    check_spectrum_integrality(A)

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
Confirm that the rotation matrix by `π/2` radians counterclockwise is not spectrum integral
(rather, it has eigenvalues `±i` [Joy15; p. 1](@cite)):
```jldoctest; setup = :(using SDiagonalizability)
julia> R = Int8.([0 -1; 1 0])
2×2 Matrix{Int8}:
 0  -1
 1   0

julia> res = check_spectrum_integrality(R);

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
eigenvalues and multiplicities of `{3: 1, -2: 4, 1: 5}` [Fox09; p. 2](@cite):
```jldoctest; setup = :(using SDiagonalizability, Graphs)
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

julia> res = check_spectrum_integrality(A);

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
If an undirected graph with integer edge weights is `{-1,0,1}`-diagonalizable (or, more
restrictively, `{-1,1}`-diagonalizable), then its Laplacian matrix has integer eigenvalues
[JP25; p. 312](@cite). Hence, validating Laplacian integrality serves as a useful screening
step in this package's principal *S*-bandwidth minimization algorithm.

It is, perhaps, an odd choice to sort the eigenvalue/multiplicity pairs in this file rather
than in [`s_bandwidth`](@ref)—after all, in the context of the overarching *S*-bandwidth
algorithm, this ordering is only ever used to determine which eigenspaces are searched for
*S*-bases first. However, the inclusion of `check_spectrum_integrality` in the public API
motivates a consistent, natural ordering of the `multiplicites` map in each (potentially
user-facing) [`SpectrumIntegralResult`](@ref) instance.

Of course, we make it a point to still sort `multiplicities` as desired (first by ascending
multiplicity then by ascending eigenvalue) in `s_bandwidth` itself; this seems more robust
than relying on `check_spectrum_integrality` to do so or even folding the logic into the
[`SpectrumIntegralResult`](@ref) inner constructor.
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
    _extract_independent_cols(A)

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
- `::AbstractMatrix{T}`: a spanning set of independent columns of `A`.

# Examples
Observe how columns with greater Euclidean norms are given priority in the pivot ordering:
```jldoctest; setup = :(using SDiagonalizability)
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
dealing with all `{-1,0,1}`-eigenvectors of a Laplacian matrix, which is the intended use
case of this helper function in this package). Our replacement tolerance, on the other hand,
is a widely accepted standard in numerical analysis which uses the maximum dimension instead
[PTVF07; p. 795](@cite).
"""
function _extract_independent_cols(A::AbstractMatrix{<:Integer})
    F = qr(A, ColumnNorm())
    rtol = _rank_rtol(A) # Use a higher tolerance (NumPy's/MATLAB's) than Julia's default

    #= In Julia 1.12+, `LinearAlgebra.rank` dispatches to a method that re-uses an existing
    QR decomposition. For compatibility with v1.10–1.11, we manually define it ourselves in
    `src/utils.jl`. =#
    r = rank(F; rtol=rtol)
    pivots = F.p[1:r] # The first `rank(A)` pivots correspond to independent columns of `A`

    return A[:, pivots]
end

"""
    _find_indices_1neg(eigvecs_01neg)

Find the indices of `{-1,1}`-eigenvectors given a map of eigenspaces.

More precisely, given a map from each eigenvalue to all `{-1,0,1}`-vectors in the associated
eigenspace, this function returns a map from each eigenvalue to the indices of the
`{-1,1}`-vectors in each value from the input map.

It is, of course, assumed that the input map indeed contains only `{-1,0,1}`-matrices as
values, and that it contains `0` as a key.

# Arguments
- `eigvecs_01neg::AbstractDict{Int,<:AbstractMatrix{Int}}`: a map from each eigenvalue to
    all `{-1,0,1}`-eigenvectors in the associated eigenspace.

# Returns
- `::AbstractDict{Int,Vector{Int}}`: a map from each eigenvalue `eigval` to the indices of
    the `{-1,1}`-columns in `eigvecs_01neg[eigval]`.

# Examples
```jldoctest; setup = :(using SDiagonalizability)
julia> eigvecs_01neg = Dict(
           0 => [1  1   1;
                 1  1   1;
                 1  1   1;
                 1  0  -1],
           2 => [ 1   0;
                  0   1;
                 -1  -1;
                  0   0],
       )
Dict{Int64, Matrix{Int64}} with 2 entries:
  0 => [1 1 1; 1 1 1; 1 1 1; 1 0 -1]
  2 => [1 0; 0 1; -1 -1; 0 0]

julia> SDiagonalizability._find_indices_1neg(eigvecs_01neg)
Dict{Int64, Vector{Int64}} with 2 entries:
  0 => [1, 3]
  2 => []
```

# Notes
It is important to note that this function implicitly assumes that the input map already
constitutes valid `{-1,0,1}`-spectra for some undirected Laplacian matrix. For instance, the
assumption that all eigenvectors corresponding to nonzero eigenvalues are orthogonal to the
all-ones vector (which must be in the kernel) is integral to the correctness of this
algorithm. Since `_find_indices_1neg` is ultimately but a helper function for
[`laplacian_spectra_01neg`](@ref), which already validates that the input passed is valid,
we choose not to double-check it here for performance reasons. (Usually, it would be good
practice to do so, but the potential sheer size of `eigvecs_01neg`, combined with the
time-sensitive nature of the algorithm, precludes the possibility in this particular case.)

Another design choice worthy of discussion is to encode (per eigenspace) the data on which
columns have entries exclusively from `{-1,1}` as a `Vector{Int}` rather than a `BitVector`.
`Vector{Int}`'s become more memory efficient than `BitVector`'s to store indices when the
density of `true`s is less than `1 / 64`, which is far greater than the proportion of
`{-1,1}`-columns in most (albeit not all) cases, simply based on our numerical results.
"""
function _find_indices_1neg(eigvecs_01neg::AbstractDict{Int,<:AbstractMatrix{Int}})
    #= The order of the Laplacian matrix. Accessing the zero key is safe, as every
    undirected graph has a rank-deficient Laplacian. =#
    n = size(eigvecs_01neg[0], 1)

    # Only even-ordered Laplacians have non-kernel {-1,1}-eigenvectors
    if n % 2 == 0 # Check all eigenspaces for {-1,1}-eigenvectors
        indices_1neg = Dict(
            eigval => findall(v -> !any(iszero, v), eachcol(eigvecs)) for
            (eigval, eigvecs) in eigvecs_01neg
        )
    else # Only check the kernel for {-1,1}-eigenvectors
        indices_1neg = Dict(eigval => Int[] for eigval in keys(eigvecs_01neg))
        indices_1neg[0] = findall(v -> !any(iszero, v), eachcol(eigvecs_01neg[0]))
    end

    return indices_1neg # Indices of the {-1,0,1}-eigenvectors without 0's
end

"""
    laplacian_spectra_01neg(L::AbstractMatrix{<:Integer})

Compute data on the `{-1,0,1}`- and `{-1,1}`-spectrum of some Laplacian matrix `L`.

The original matrix, whether it is `{-1,0,1}`- and `{-1,1}`-diagonalizable, and further data
on the eigenvalues and their corresponding `{-1,0,1}`- and `{-1,1}`-eigenvectors are all
stored in a [`LaplacianSpectrum01Neg`](@ref) instance and returned.

# Arguments
- `L::AbstractMatrix{<:Integer}`: the Laplacian matrix in whose `{-1,0,1}`- and
    `{-1,1}`-spectrum we are interested.

# Returns
- `::LaplacianSpectrum01Neg`: a struct containing the following fields:
    - `laplacian_matrix`: the `L` matrix, copied to avoid shared mutability.
    - `diagonalizable_01neg`: whether the Laplacian is `{-1,0,1}`-diagonalizable.
    - `diagonalizable_1neg`: whether the Laplacian is `{-1,1}`-diagonalizable.
    - `multiplicities`: a map from each eigenvalue to its multiplicity, sorted first by
        ascending multiplicity then by ascending eigenvalue. (This field is `nothing` if the
        Laplacian is not spectrum integral.)
    - `eigvecs_01neg`: a map from each eigenvalue to its `{-1,0,1}`-eigenvectors. (This
        field is `nothing` if the Laplacian is not spectrum integral.)
    - `indices_1neg`: a map from each eigenvalue to the indices of the `{-1,1}`-columns in
        `eigvecs_01neg[eigval]`. (This field is `nothing` if the Laplacian is not spectrum
        integral.)
    - `bases_01neg`: a map from each eigenvalue to a `{-1,0,1}`-basis for the corresponding
        eigenspace. (This field is `nothing` if the Laplacian is not spectrum integral.)
    - `bases_1neg`: a map from each eigenvalue to a `{-1,1}`-basis for the corresponding
        eigenspace. (This field is `nothing` if the Laplacian is not spectrum integral.)

# Examples
TODO: Add examples

# Notes
`laplacian_spectra_01neg` outsources the actual computation logic to the
[`_typed_laplacian_spectra_01neg`](@ref) helper, which dispatches to separate methods for
different categories of Laplacians. The input matrix is first classified and cast to a
concrete subtype of [`_TypedLaplacian`](@ref), allowing said helper to apply specialized
algorithms relying on the different `{-1,0,1}`-spectral properties of different types of
graphs. See the [`_TypedLaplacian`](@ref) documentation for more details on the different
classifications of Laplacians and the [`_typed_laplacian_spectra_01neg`](@ref) documentation
for more details on the optimized algorithm used for each type.
"""
function laplacian_spectra_01neg(L::AbstractMatrix{<:Integer})
    TL = _cast_to_typed_laplacian(L)
    return _typed_laplacian_spectra_01neg(TL)
end

"""
    _typed_laplacian_spectra_01neg(TL::_TypedLaplacian)

TODO: Write here
"""
function _typed_laplacian_spectra_01neg(TL::_TypedLaplacian)
    throw(NotImplementedError(_typed_laplacian_spectra_01neg, typeof(TL), _TypedLaplacian))
end

"""
    _typed_laplacian_spectra_01neg(TL::_NullGraphLaplacian)

TODO: Write here
"""
function _typed_laplacian_spectra_01neg(TL::_NullGraphLaplacian)
    eigspaces_01neg = OrderedDict{Int,_Eigenspace01Neg}()
    return LaplacianSpectrum01Neg(TL.matrix, true, true, eigspaces_01neg)
end

"""
    _typed_laplacian_spectra_01neg(TL::_EmptyGraphLaplacian)

TODO: Write here
"""
function _typed_laplacian_spectra_01neg(TL::_EmptyGraphLaplacian)
    L = TL.matrix
    n = size(L, 1)

    #= JET.jl cannot infer whether `kernel_eigvecs_01neg` is a `Matrix{Int}` or a
    `Vector{Any}`, so we explicitly assert its type here for compiler inference during
    static analysis. =#
    kernel_eigvecs_01neg = hcat(_pot_kernel_eigvecs_01neg(n)...)
    @assert kernel_eigvecs_01neg isa Matrix{Int}

    eigvecs_01neg = Dict(0 => kernel_eigvecs_01neg)
    indices_1neg = _find_indices_1neg(eigvecs_01neg)
    pot_kernel_basis_1neg = _extract_independent_cols(kernel_eigvecs_01neg)

    if size(pot_kernel_basis_1neg, 2) == n # The potential {-1,1}-basis spans the kernel
        kernel_basis_01neg = pot_kernel_basis_1neg
        kernel_basis_1neg = copy(pot_kernel_basis_1neg) # Avoid shared mutability
        diagonalizable_1neg = true
    else # No {-1,1}-basis exists, so compute a {-1,0,1}-eigenbasis (guaranteed to exist)
        kernel_basis_01neg = _extract_independent_cols(kernel_eigvecs_01neg)
        kernel_basis_1neg = nothing
        diagonalizable_1neg = false
    end

    eigspaces_01neg = OrderedDict(
        0 => _Eigenspace01Neg(
            n,
            kernel_eigvecs_01neg,
            indices_1neg[0],
            kernel_basis_01neg,
            kernel_basis_1neg,
        ),
    )
    return LaplacianSpectrum01Neg(L, true, diagonalizable_1neg, eigspaces_01neg)
end

"""
    _typed_laplacian_spectra_01neg(TL::_CompleteGraphLaplacian)

TODO: Write here
"""
function _typed_laplacian_spectra_01neg(TL::_CompleteGraphLaplacian)
    L = TL.matrix
    n = size(L, 1)

    #= For every `n ≥ 2`, `Kₙ` (assuming a uniform edge weight of `w`) has eigenvalues 0
    with multiplicity 1 and `n * w` with multiplicity `n - 1`. (`n` is guaranteed to be at
    least 2 in this context due to the Laplacian type casting logic.) =#
    eigval_nonzero = n * TL.weight

    kernel_eigvecs_01neg = ones(Int, n, 1) # The all-ones vector spans the kernel

    #= The non-kernel eigenspace necessarily contains all vectors orthogonal to the kernel,
    so we do not need to filter through the output of the generator. =#
    #= Additionally, JET.jl cannot infer whether `nonkernel_eigvecs_01neg` is a
    `Matrix{Int}` or a `Vector{Any}`, so we explicitly assert its type here for compiler
    inference during static analysis. =#
    nonkernel_eigvecs_01neg = hcat(_pot_nonkernel_eigvecs_01neg(n)...)
    @assert nonkernel_eigvecs_01neg isa Matrix{Int}

    # All the {-1,0,1}-eigenvectors corresponding to each eigenvalue, as columns of a matrix
    eigvecs_01neg = Dict(
        0 => kernel_eigvecs_01neg, eigval_nonzero => nonkernel_eigvecs_01neg
    )
    indices_1neg = _find_indices_1neg(eigvecs_01neg)

    # Again, the set containing the all-ones vector forms a basis for the kernel
    kernel_basis_01neg = copy(kernel_eigvecs_01neg) # Avoid shared mutability
    kernel_basis_1neg = copy(kernel_eigvecs_01neg) # Avoid shared mutability
    pot_nonkernel_basis_1neg = _extract_independent_cols(
        nonkernel_eigvecs_01neg[:, indices_1neg[eigval_nonzero]]
    )

    #= Find a {-1,0,1}- and, if one exists, a {-1,1}-basis for the remaining eigenspace.
    {-1,1}-bases are preferable, so they take priority whenever they exist. =#
    if size(pot_nonkernel_basis_1neg, 2) == n - 1 # The potential {-1,1}-basis spans the space
        nonkernel_basis_01neg = pot_nonkernel_basis_1neg
        nonkernel_basis_1neg = copy(pot_nonkernel_basis_1neg) # Avoid shared mutability
        diagonalizable_1neg = true
    else # No {-1,1}-basis exists, so use a {-1,0,1}-eigenbasis (guaranteed to exist). =#
        nonkernel_basis_01neg = _extract_independent_cols(nonkernel_eigvecs_01neg)
        nonkernel_basis_1neg = nothing
        diagonalizable_1neg = false
    end

    eigspaces_01neg = OrderedDict{Int,_Eigenspace01Neg}()
    kernel_01neg = _Eigenspace01Neg(
        n, kernel_eigvecs_01neg, indices_1neg[0], kernel_basis_01neg, kernel_basis_1neg
    )
    nonkernel_01neg = _Eigenspace01Neg(
        n - 1,
        nonkernel_eigvecs_01neg,
        indices_1neg[eigval_nonzero],
        nonkernel_basis_01neg,
        nonkernel_basis_1neg,
    )

    # Sort first by ascending multiplicity then by ascending eigenvalue
    if n == 2 && eigval_nonzero < 0
        eigspaces_01neg[eigval_nonzero] = nonkernel_01neg
        eigspaces_01neg[0] = kernel_01neg
    else
        eigspaces_01neg[0] = kernel_01neg
        eigspaces_01neg[eigval_nonzero] = nonkernel_01neg
    end

    return LaplacianSpectrum01Neg(L, true, diagonalizable_1neg, eigspaces_01neg)
end

"""
    _typed_laplacian_spectra_01neg(TL::_ArbitraryGraphLaplacian)

TODO: Write here
"""
function _typed_laplacian_spectra_01neg(TL::_ArbitraryGraphLaplacian)
    L = TL.matrix
    res = check_spectrum_integrality(L)

    if !res.spectrum_integral
        # For integer edge weights, {-1,0,1}-diagonalizability implies Laplacian integrality
        return LaplacianSpectrum01Neg(L, false, false, nothing)
    end

    #= `multiplicities` is only `nothing` when `spectrum_integral` is false, so we
    explicitly assert its type here for compiler inference during static analysis. =#
    multiplicities = res.multiplicities
    @assert !isnothing(multiplicities)

    n = size(L, 1)
    eigvals_nonzero = filter(!iszero, keys(multiplicities))
    eigvecs_01neg = Dict{Int,AbstractMatrix{Int}}(
        # Initialize resizeable matrices to store an unknown number of eigenvectors
        eigval => ElasticMatrix{Int}(undef, n, 0) for eigval in eigvals_nonzero
    )

    # The kernel is simpler and handled differently from the other eigenspaces
    if multiplicities[0] == 1 # The all-ones vector spans the kernel
        eigvecs_01neg[0] = ones(Int, n, 1)
    else # Check all {-1,0,1}-vectors in ℝⁿ, unique up to span
        # TODO: Multithread
        eigvecs_01neg[0] = hcat(
            Iterators.filter(v -> iszero(L * v), _pot_kernel_eigvecs_01neg(n))...
        )
    end

    # TODO: Multithread
    # Now fill up the remaining (non-kernel) eigenspaces
    for v in _pot_nonkernel_eigvecs_01neg(n)
        for eigval in eigvals_nonzero
            if L * v == eigval * v
                append!(eigvecs_01neg[eigval], v)
                break
            end
        end
    end

    indices_1neg = _find_indices_1neg(eigvecs_01neg)

    # TODO: Comment on why we type the values
    bases_01neg = Dict{Int,Union{Nothing,Matrix{Int}}}(
        eigval => _extract_independent_cols(vecs) for (eigval, vecs) in eigvecs_01neg
    )
    bases_1neg = Dict{Int,Union{Nothing,Matrix{Int}}}(
        eigval => _extract_independent_cols(vecs[:, indices_1neg[eigval]]) for
        (eigval, vecs) in eigvecs_01neg
    )
    diagonalizable_01neg = true
    diagonalizable_1neg = true

    #= If they exist, find a {-1,0,1}- and {-1,1}-basis for each non-kernel eigenspace.
    {-1,1}-bases are preferable, so they take priority whenever they exist. =#
    for eigval in eigvals_nonzero
        # TODO: Comment on the `@assert`'s for compiler inference during static analysis
        pot_basis_01neg = bases_01neg[eigval]
        @assert !isnothing(pot_basis_01neg)
        pot_basis_1neg = bases_1neg[eigval]
        @assert !isnothing(pot_basis_1neg)

        multiplicity = multiplicities[eigval]

        if size(pot_basis_01neg, 2) < multiplicity # The eigenspace has no {-1,0,1}-basis
            bases_01neg[eigval] = bases_1neg[eigval] = nothing
            diagonalizable_01neg = diagonalizable_1neg = false
            # There exists a {-1,0,1}- but not a {-1,1}-eigenbasis
        elseif size(pot_basis_1neg, 2) < multiplicity
            bases_1neg[eigval] = nothing
            diagonalizable_1neg = false
        else # There exists a {-1,1}-eigenbasis, so use it in both dictionaries
            bases_01neg[eigval] = copy(pot_basis_1neg) # Avoid shared mutability
        end
    end

    # Create and populate our usual wrapper for the eigenspace data
    eigspaces_01neg = OrderedDict{Int,_Eigenspace01Neg}()

    for (eigval, multiplicity) in multiplicities
        eigspaces_01neg[eigval] = _Eigenspace01Neg(
            multiplicity,
            eigvecs_01neg[eigval],
            indices_1neg[eigval],
            bases_01neg[eigval],
            bases_1neg[eigval],
        )
    end

    return LaplacianSpectrum01Neg(
        L, res.spectrum_integral, diagonalizable_1neg, eigspaces_01neg
    )
end
