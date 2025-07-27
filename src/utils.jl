# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    _assert_matrix_is_undirected_laplacian(L) -> Nothing

Validate that `L` is the Laplacian matrix of an undirected, possibly weighted graph.

# Arguments
- `L::AbstractMatrix{<:Integer}`: the purported Laplacian matrix.

# Returns
- `nothing`: if the check is passed, no output is produced.

# Throws
- `DomainError`: if `L` is not symmetric or has nonzero row sums.

# Examples
The Laplacian matrix of the (undirected) star graph on ``5`` vertices passes the check:
```jldoctest
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

The adjacency matrix of the (undirected) cycle graph on ``4`` vertices is symmetric but has
nonzero row sums, so it fails the check:
```jldoctest
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
```jldoctest
julia> using Graphs

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
sums, since the ``(i, i)``-th entry is the weighted degree of node ``i`` (the sum of all
incident edges' weights) and the ``(i, j)``-th entry for ``i ≠ j`` is the negation of the
weight of edge ``(i, j)`` (or simply ``0``, if no such edge exists).

Given the highly optimized, lazy, zero-allocation implementation of
`LinearAlgebra.issymmetric`, the symmetry check is performed first. (Both steps are
``O(n²)`` in the worst case, but testing for symmetry is far more performant in practice.)
This also allows us to (also lazily) check for nonzero column sums rather than nonzero row
sums (since these are equivalent for symmetric matrices) in the second step, taking
advantage of Julia's column-major storage model.

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
    function _rank_rtol(A::AbstractMatrix{<:Real}) -> Float64

Return a reasonable relative tolerance for computing matrix rank via SVD or QRD.

The output is intended to be passed as a keyword argument to `LinearAlgebra.rank`.
`LinearAlgebra.eigtype` is used to determine the appropriate machine epsilon for the element
type of `A` when `eltype(A)` is not an `AbstractFloat`.

# Arguments
- `A::AbstractMatrix{<:Real}`: the matrix for which to compute a tolerance.

# Returns
- `tol::Float64`: a reasonable relative tolerance for computing matrix rank via SVD or QRD.
    This scales proportionally to the maximum dimension of `A`.

# Notes
`LinearAlgebra.rank`'s default `rtol` of `min(m,n) * ϵ` for computing the rank of an
``m×n`` matrix may result in overestimating rank when ``|m - n| ≫ 0``, since condition
number (which determines how numerically stable SVD and QRD are) grows with both dimensions
[CD05; p. 603](@cite). Given that we often deal with short-and-fat matrices in this package
(particularly when processing all ``{-1,0,1}``-eigenvectors of a Laplacian matrix), we turn
instead to the same relative tolerance used by NumPy's and MATLAB's rank
functions—`max(m,n) * ϵ` [Num25, MAT25](@cite). (Indeed, this is a widely adopted standard
across the field of numerical analysis [PTVF07; p. 795](@cite).)
"""
function _rank_rtol(A::AbstractMatrix{T}) where {T}
    return maximum(size(A)) * eps(LinearAlgebra.eigtype(T))
end

#= In Julia 1.12+, `LinearAlgebra.rank` dispatches to a method that re-uses an existing QR
decomposition. For compatibility with v1.10–1.11, we manually define it ourselves here. =#
@static if VERSION < v"1.12"
    #! format: off
    function LinearAlgebra.rank(A::QRPivoted; atol::Real=0, rtol::Real=min(size(A)...) * eps(real(float(eltype(A)))) * iszero(atol))
        m = min(size(A)...)
        m == 0 && return 0
        tol = max(atol, rtol*abs(A.factors[1,1]))
        return something(findfirst(i -> abs(A.factors[i,i]) <= tol, 1:m), m+1) - 1
    end
    #! format: on

    @doc """
        rank(A::QRPivoted{<:Any, T}; atol::Real=0, rtol::Real=min(n,m)*ϵ) where {T}

    Compute the numerical rank of the QR factorization `A` by counting how many diagonal entries of
    `A.factors` are greater than `max(atol, rtol*Δ₁)` where `Δ₁` is the largest calculated such entry.
    This is similar to the `rank(::AbstractMatrix)` method insofar as it counts the number of
    (numerically) nonzero coefficients from a matrix factorization, although the default method uses an
    SVD instead of a QR factorization. Like `rank(::SVD)`, this method also re-uses an existing
    matrix factorization.

    Computing rank via QR factorization should almost always produce the same results as via SVD,
    although this method may be more prone to overestimating the rank in pathological cases where the
    matrix is ill-conditioned. It is also worth noting that it is generally faster to compute a QR
    factorization than an SVD, so this method may be preferred when performance is a concern.

    `atol` and `rtol` are the absolute and relative tolerances, respectively.
    The default relative tolerance is `n*ϵ`, where `n` is the size of the smallest dimension of `A`
    and `ϵ` is the `eps` of the element type of `A`.

    !!! note
        When accessed directly via `LinearAlgebra`, the `rank(::QRPivoted)` method requires at least
        Julia 1.12, so `SDiagonalizability` defines this method manually for compatibility with
        v1.10–1.11.
    """ -> rank
end
