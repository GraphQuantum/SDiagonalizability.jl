# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    _sort_tuple(S::T) where {T<:Tuple{Vararg{Integer}}} -> T

Return a sorted version of the input tuple.

# Arguments
- `S::T<:Tuple{Vararg{Integer}}`: the tuple to sort.

# Returns
- `::T`: a sorted version of `S`.

# Examples
```jldoctest
julia> S = (-1, 1, 0)
(-1, 1, 0)

julia> SDiagonalizability._sort_tuple(S)
(-1, 0, 1)
```
"""
function _sort_tuple(S::Tuple{Vararg{Integer}})
    return Tuple(sort(collect(S)))
end

"""
    _extract_independent_cols(A) -> Matrix{Int}

Return a (not necessarily unique) independent spanning subset of the columns of `A`.

Computing a rank-revealing (pivoted) QR decomposition of `A`, the scaling coefficients from
the orthogonalization process are used to determine the rank (rather than recompute it with
an SVD), while the pivots are used to extract a spanning set of independent columns.

The rank-revealing Businger–Golub QR algorithm is used for the pivoting strategy, appending
the "most independent" column with respect to the current set of pivots at each step via
Householder transformations [BG65; pp. 269--70].

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
compute the rank of `A` rather than compute a separate SVD. We thus count the nonzero
scaling coefficients—that is, the diagonal entries of the `R` matrix in `A = QR`—to
determine the rank, similarly to how we count the nonzero singular values in an SVD.

It is worth noting that we manually specify a higher relative tolerance for this rank
computation. Further discussion can be found in the [`_rank_rtol`](@ref) documentation, but
in short, a critical part of the formula for `LinearAlgebra.rank`'s default `rtol`
uses the minimum dimension of the input matrix. This may result in rank overestimation for
tall-and-skinny and short-and-fat matrices (precisely the type we expect to encounter when
dealing with all ``\\{-1, 0, 1\\}``-eigenvectors of a Laplacian matrix, which is the
intended use case of this helper function in this package). Our replacement tolerance, on
the other hand, is a widely accepted standard in numerical analysis which uses the maximum
dimension instead [PTVF07; p. 795].

# References

- [BG65](@cite): P. Businger and G. H. Golub. *Linear Least Squares Solutions by Householder
    Transformations*. *Numerische Mathematik* **7**, 269–76 (1965).
    https://doi.org/10.1007/BF01436084.

- [PTVF07](@cite): W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery.
    *Numerical Recipes: The Art of Scientific Computing*. 3rd Edition (Cambridge University
    Press, Cambridge, UK, 2007). ISBN: 978-0-521-88068-8.
    https://dl.acm.org/doi/10.5555/1403886.
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
ways of extending the concept of the Laplacian to directed graphs [VL20; p. 196].)
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

# References

- [VL20](@cite): J. J. Veerman and R. Lyons. *A Primer on Laplacian Dynamics in Directed
    Graphs*. *Nonlinear Phenomena in Complex Systems* **23**, 196–206 (2020).
    https://doi.org/10.33581/1561-4085-2020-23-2-196-206.
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
    _assert_graph_has_defined_s_bandwidth(g) -> Nothing

[TODO: Write here]
"""
function _assert_graph_has_defined_s_bandwidth(g::AbstractGraph)
    if is_directed(g)
        throw(DomainError(g, "*S*-bandwidth is not defined for directed graphs"))
    end

    if has_self_loops(g)
        throw(
            DomainError(
                graph,
                "*S*-bandwidth is not defined for multigraphs; got a graph with self-loops",
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
[CD05; p. 603]. Given that we often deal with short-and-fat matrices in this package
(particularly when processing all ``\\{-1, 0, 1\\}``-eigenvectors of a Laplacian matrix), we
turn instead to the same relative tolerance used by NumPy's and MATLAB's rank
functions—`max(m,n) * ϵ` [Num25, MAT25]. (Indeed, this is a widely adopted standard across
the field of numerical analysis [PTVF07; p. 795].)

# References

- [CD05](@cite): Z. Chen and J. Dongarra. *Condition Numbers of Gaussian Random Matrices*.
    *SIAM Journal on Matrix Analysis and Applications* **27**, 603–20 (2005).
    https://doi.org/10.1137/040616413.
- [MAT25](@cite): MATLAB Developers, *rank*. MATLAB reference documentation – R2025a (2025).
    Accessed: 2025-05-29. https://www.mathworks.com/help/matlab/ref/rank.html.
- [Num25](@cite): NumPy Developers, *numpy.linalg.matrix_rank*. NumPy reference
    documentation – v2.2 (2025). Accessed: 2025-05-22.
    https://numpy.org/doc/stable/reference/generated/numpy.linalg.matrix_rank.html.
- [PTVF07](@cite): W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery.
    *Numerical Recipes: The Art of Scientific Computing*. 3rd Edition (Cambridge University
    Press, Cambridge, UK, 2007). ISBN: 978-0-521-88068-8.
    https://dl.acm.org/doi/10.5555/1403886.
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
