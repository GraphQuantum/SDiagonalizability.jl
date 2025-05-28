# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

#= TODO: Extend functionality to multi-argument functions. Talk about "dispatch to [type] for
# [parameter]" rather than the current error message. =#
struct NotImplementedError <: Exception
    f::Function
    concretetype::Type
    abstracttype::Type

    function NotImplementedError(f::Function, concretetype::Type, abstracttype::Type)
        if !isconcretetype(concretetype)
            throw(ArgumentError("Expected a concrete type, got $concretetype"))
        end

        if !isabstracttype(abstracttype)
            throw(ArgumentError("Expected an abstract type, got $abstracttype"))
        end

        if !(concretetype <: abstracttype)
            throw(ArgumentError("Expected a subtype of $abstracttype, got $concretetype"))
        end

        return new(f, concretetype, abstracttype)
    end
end

function Base.showerror(io::IO, e::NotImplementedError)
    print(
        io,
        """NotImplementedError with `$(e.concretetype)`:
        The function `$(e.f)` is not yet implemented for this subtype of `$(e.abstracttype)`.""",
    )
end

"""
    _assert_matrix_is_undirected_laplacian(L)

Validate that `L` is the Laplacian matrix of an undirected, possibly weighted graph.

# Arguments
- `L::AbstractMatrix{<:Integer}`: the purported Laplacian matrix.

# Returns
- `nothing`: if the check is passed, no output is produced.

# Throws
- `DomainError`: if `L` is not symmetric or has nonzero row sums.

# Examples
The Laplacian matrix of the (undirected) star graph on `5` vertices passes the check:
```jldoctest; setup = :(using SDiagonalizability)
julia> L = [ 4  -1  -1  -1  -1;
            -1   1   0   0   0;
            -1   0   1   0   0;
            -1   0   0   1   0;
            -1   0   0   0   1]
5×5 Matrix{Int64}:
  4  -1  -1  -1  -1
 -1   1   0   0   0
 -1   0   1   0   0
 -1   0   0   1   0
 -1   0   0   0   1

julia> isnothing(SDiagonalizability._assert_matrix_is_undirected_laplacian(L))
true
```

The adjacency matrix of the (undirected) cycle graph on `4` vertices is symmetric but has
nonzero row sums, so it fails the check:
```jldoctest; setup = :(using SDiagonalizability)
julia> A = [0  1  0  1;
            1  0  1  0;
            0  1  0  1;
            1  0  1  0]
4×4 Matrix{Int64}:
 0  1  0  1
 1  0  1  0
 0  1  0  1
 1  0  1  0

julia> SDiagonalizability._assert_matrix_is_undirected_laplacian(A)
ERROR: DomainError with [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0]:
Matrix has nonzero row sums; cannot be an (undirected) Laplacian
[...]
```

Both the in-degree and out-degree Laplacian matrices of this random tournament digraph have
zero row sums but are not symmetric, so they fail the check. (These are the two standard
ways of extending the concept of the Laplacian to directed graphs [VL20; p. 196](@cite).)
```jldoctest; setup = :(using SDiagonalizability, Graphs)
julia> G = random_tournament_digraph(3; seed=87)
{3, 3} directed simple Int64 graph

julia> L_in = laplacian_matrix(G; dir=:in)
3×3 SparseArrays.SparseMatrixCSC{Int64, Int64} with 5 stored entries:
  ⋅  ⋅   ⋅
 -1  2  -1
 -1  ⋅   1

julia> SDiagonalizability._assert_matrix_is_undirected_laplacian(L_in)
ERROR: DomainError with sparse([2, 3, 2, 2, 3], [1, 1, 2, 3, 3], [-1, -1, 2, -1, 1], 3, 3):
Matrix is not symmetric; cannot be an (undirected) Laplacian
[...]

julia> L_out = laplacian_matrix(G; dir=:out)
3×3 SparseArrays.SparseMatrixCSC{Int64, Int64} with 5 stored entries:
 2  -1  -1
 ⋅   ⋅   ⋅
 ⋅  -1   1

julia> SDiagonalizability._assert_matrix_is_undirected_laplacian(L_out)
ERROR: DomainError with sparse([1, 1, 3, 1, 3], [1, 2, 2, 3, 3], [2, -1, -1, -1, 1], 3, 3):
Matrix is not symmetric; cannot be an (undirected) Laplacian
[...]
```

# Notes
If edges are to be bidirectional, then `L` must be symmetric. `L` must also have zero row
sums, since the `(i, i)`-th entry is the weighted degree of node `i` (the sum of all
incident edges' weights) and the `(i, j)`-th entry for `i ≠ j` is the negation of the weight
of edge `(i, j)` (or simply `0`, if no such edge exists).

Given the highly optimized, lazy, zero-allocation implementation of
`LinearAlgebra.issymmetric`, the symmetry check is performed first. (Both steps are `O(n²)`
in the worst case, but testing for symmetry is far more performant in practice.) This also
allows us to (also lazily) check for nonzero column sums rather than nonzero row sums (since
these are equivalent for symmetric matrices) in the second step, taking advantage of Julia's
column-major storage model.

At first blush, it may seem as though the choice of `DomainError` over something like
`ArgumentError` (or even simply the return of a boolean) constitutes poor design. However,
this is informed by the simple *ad hoc* use of this function to validate inputs for other
functions requiring Laplacian matrices. Certainly, this function is never meant to be
publicly exposed on its own.
"""
function _assert_matrix_is_undirected_laplacian(L::AbstractMatrix{<:Integer})
    if !issymmetric(L)
        throw(
            DomainError(L, "Matrix is not symmetric; cannot be an (undirected) Laplacian")
        )
    end

    # This method of lazy evaluation is faster than calling `!iszero(sum(L, dims=1))`
    if any(col -> !iszero(sum(col)), eachcol(L))
        throw(
            DomainError(
                L, "Matrix has nonzero row sums; cannot be an (undirected) Laplacian"
            ),
        )
    end

    return nothing
end

"""
    function _rank_rtol(A::AbstractMatrix{<:Real})

Return a reasonable relative tolerance for computing matrix rank via SVD or QRD.

If using SVD to compute rank with singular values Ω, then the intended use case is to
numerically estimate the rank `r` of `A` by `r = count(Ω .> maximum(Ω) * _rank_rtol(A))`.
Similarly, if using QRD to compute rank with scaling coefficients Δ (the diagonal entries of
the `R` matrix in `A = QR`), then the user should numerically estimate the rank by `r =
count(Δ .> maximum(Δ) * _rank_rtol(A))`.

Since the machine epsilon is only defined for floating-point types, we call `eps()`
(defaulting to `Float64`) instead of `eps(eltype(A))` when computing the tolerance. If the
elements of `A` are floats, on the other hand, then this function dispatches to a separate
method that uses the element type's machine epsilon.

# Arguments
- `A::AbstractMatrix{<:Real}`: the matrix for which to compute a tolerance.

# Returns
- `tol::Float64`: a reasonable relative tolerance for computing matrix rank via SVD or QRD.
    This scales proportionally to the maximum dimension of `A`.

# Notes
This is very much an *ad hoc* function applied specifically for numerical stability within
this package and is never meant to be used as part of the public API. The
`sqrt(maximum(size(A)) * eps())` formula is inspired by NumPy's default relative tolerance
of `maximum(size(A)) * eps()` in the `numpy.linalg.matrix_rank` function [numpy25](@cite),
with a square root taken for additional robustness. (In particular, the short, fat matrices
holding all `{-1,0,1}`-eigenvectors of a Laplacian matrix are often ill-conditioned, hence
the need for this higher tolerance.)

In the case of using SVD to compute matrix rank, the decision to default to the `Float64`
machine epsilon is motivated by how `LinearAlgebra.rank` handles non-floating-point
matrices, as LAPACK automatically converts to `Float64` under the hood. In the case of QRD,
our *ad hoc* use of this function in [`_extract_independent_cols`](@ref) takes in
floating-point matrices anyway, so this particular method is no longer relevant.
"""
function _rank_rtol(A::AbstractMatrix{<:Real})
    return sqrt(maximum(size(A)) * eps())
end

"""
    function _rank_rtol(A::AbstractMatrix{<:AbstractFloat})

Return a reasonable relative tolerance for computing matrix rank via SVD or QRD.

If using SVD to compute rank with singular values `Ω`, then the intended use case is to
numerically estimate the rank `r` of `A` by `r = count(Ω .> maximum(Ω) * _rank_rtol(A))`.
Similarly, if using QRD to compute rank with scaling coefficients `Δ` (the diagonal entries
of the `R` matrix in `A = QR`), then the user should numerically estimate the rank by `r =
count(Δ .> maximum(Δ) * _rank_rtol(A))`.

# Arguments
- `A::AbstractMatrix{T<:AbstractFloat}`: the matrix for which to compute a tolerance.

# Returns
- `tol::Float64`: a reasonable relative tolerance for computing matrix rank via SVD or QRD.
    This scales proportionally to the maximum dimension of `A` and the machine epsilon of
    the element type `T`.

# Notes
This is very much an *ad hoc* function applied specifically for numerical stability within
this package and is never meant to be used as part of the public API. The
`sqrt(maximum(size(A)) * eps())` formula is inspired by NumPy's default relative tolerance
of `maximum(size(A)) * eps()` in the `numpy.linalg.matrix_rank` function [numpy25](@cite),
with a square root taken for additional robustness. (In particular, the short, fat matrices
holding all `{-1,0,1}`-eigenvectors of a Laplacian matrix are often ill-conditioned, hence
the need for this higher tolerance.)
"""
function _rank_rtol(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    return sqrt(maximum(size(A)) * Float64(eps(T)))
end
