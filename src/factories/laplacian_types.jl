# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    abstract type _TypedLaplacian

An abstract type representing a classified Laplacian matrix of an undirected graph.

# Interface
Concrete subtypes of `_TypedLaplacian` **must** implement the following fields:
- `matrix::Matrix{Int}`: the Laplacian matrix of the graph.

# Subtypes
- [_NullGraphLaplacian`](@ref): a Laplacian wrapper for the (unique) graph with no nodes.
- [`_EmptyGraphLaplacian`](@ref): a Laplacian wrapper for any graph with no edges on `n ≥ 1`
    nodes.
- [`_CompleteGraphLaplacian`](@ref): a Laplacian wrapper for any graph on `n ≥ 2` nodes
    where every pair of nodes is connected by an edge and all edges possess the same weight.
- [`_ArbitraryGraphLaplacian`](@ref): a Laplacian wrapper for any graph on `n ≥ 3` nodes
    with no particular strcuture of relevance in the context of investigating
    `{-1,0,1}`-spectra.

# Notes
Given that this type and all its subtypes are privately encapsulated, we expect
`_TypedLaplacian` instances to be created exclusively via the
[`_cast_to_typed_laplacian`](@ref) factory function, which copies each input matrix to avoid
shared mutability and converts it to a `Matrix{Int}` in the process. Therefore, it is safe
to restrict the type of `matrix` field (in each concrete subtype) to `Matrix{Int}` in lieu
of the more generic `AbstractMatrix{<:Integer}`. This decision may contribute marginal
performance benefits via reduced type complexity and is highly unlikely to cause any issues
(as this type is never to be exported to the public API). It also serves as an additional
safeguard against directly passing a (wrongly typed) matrix to a subtype's constructor
rather than calling [`_cast_to_typed_laplacian`](@ref).

On a related note, it seems prudent to discuss the intended use case of `_TypedLaplacian`.
The interface and its subtypes are meant to support [`_cast_to_typed_laplacian`](@ref). This
factory, in turn, determines which algorithm [`laplacian_spectra_01neg`](@ref) is to
implement when computing data on the Laplacian matrix of a given graph, as different types
of graphs exhibit different spectral properties which we may exploit to improve performance.
"""
abstract type _TypedLaplacian end

"""
    struct _NullGraphLaplacian

A wrapper for the Laplacian matrix of the (unique) graph with no nodes.

# Fields
- `matrix::Matrix{Int}`: the Laplacian in question. (Necessarily a `0×0 Matrix{Int}`.)

# Supertype Hierarchy
`_NullGraphLaplacian` <: [`_TypedLaplacian`](@ref) <: Any

# Notes
See the [`_typed_laplacian_spectra_01neg(::_NullGraphLaplacian)`](@ref) documentation for
more details on the `{-1,0,1}`-spectral properties of this Laplacian type.

See also the documentation for parent type [`_TypedLaplacian`](@ref) if seeking
justification for typing the `matrix` field as `Matrix{Int}` rather than
`AbstractMatrix{<:Integer}`.
"""
struct _NullGraphLaplacian <: _TypedLaplacian
    matrix::Matrix{Int}
end

"""
    struct _EmptyGraphLaplacian

A wrapper for the Laplacian matrix of any graph with no edges on `n ≥ 1` nodes.

# Fields
- `matrix::Matrix{Int}`: the `n`-by-`n` Laplacian in question. (Necessarily equal to
    `zeros(Int, n, n)` for some `n ≥ 1`.)

# Supertype Hierarchy
`_EmptyGraphLaplacian` <: [`_TypedLaplacian`](@ref) <: Any

# Notes
We are confident in the assumption that `n` is at least `1` because the only graph on zero
nodes is the null graph, which is handled by the [`_NullGraphLaplacian`](@ref) type.

See the [`_typed_laplacian_spectra_01neg::_EmptyGraphLaplacian`](@ref) documentation for
more details on the `{-1,0,1}`-spectral properties of this Laplacian type.

See also the documentation for parent type [`_TypedLaplacian`](@ref) if seeking
justification for typing the `matrix` field as `Matrix{Int}` rather than
`AbstractMatrix{<:Integer}`.
"""
struct _EmptyGraphLaplacian <: _TypedLaplacian
    matrix::Matrix{Int}
end

"""
    struct _CompleteGraphLaplacian

A wrapper for the Laplacian of any uniformly weighted complete graph on `n ≥ 2` nodes.

# Fields
- `matrix::Matrix{Int}`: the `n`-by-`n` Laplacian in question. (Necessarily symmetric with
    zero row sums and all off-diagonal entries equal to `-weight` for some `n ≥ 2`.)
- `weight::Int`: the weight of each edge in the graph represented by `matrix`. (Necessarily
    a nonzero integer.)

# Supertype Hierarchy
`_CompleteGraphLaplacian` <: [`_TypedLaplacian`](@ref) <: Any

# Notes
The only graph on zero nodes is the null graph (handled by the [`_NullGraphLaplacian`](@ref)
type) and the only graph on one node is an empty graph (handled by the
[`_EmptyGraphLaplacian`](@ref) type), so we are confident in the assumption that `n ≥ 2` for
any `_CompleteGraphLaplacian` instance.

See the [`_typed_laplacian_spectra_01neg::_CompleteGraphLaplacian`](@ref) documentation for
more details on the `{-1,0,1}`-spectral properties of this Laplacian type.

See also the documentation for parent type [`_TypedLaplacian`](@ref) if seeking
justification for typing the `matrix` field as `Matrix{Int}` rather than
`AbstractMatrix{<:Integer}`.
"""
struct _CompleteGraphLaplacian <: _TypedLaplacian
    matrix::Matrix{Int}
    weight::Int
end

"""
    struct _ArbitraryGraphLaplacian

A wrapper for the Laplacian of any graph on `n ≥ 3` nodes with no particular structure
of interest.

# Fields
- `matrix::Matrix{Int}`: the `n`-by-`n` Laplacian in question. (Necessarily symmetric
    with zero row sums for some `n ≥ 3`.)

# Supertype Hierarchy
`_ArbitraryGraphLaplacian` <: [`_TypedLaplacian`](@ref) <: Any

# Notes
The only graph on zero nodes is the null graph (handled by the [`_NullGraphLaplacian`](@ref)
type), the only graph on one node is an empty graph (handled by the
[`_EmptyGraphLaplacian`](@ref) type), and the only graphs on two nodes are an empty graph
and a complete graph (handled by the [`_CompleteGraphLaplacian`](@ref) type), so we are
confident in the assumption that `n ≥ 3` for any `_ArbitraryGraphLaplacian` instance.

See the [`_typed_laplacian_spectra_01neg::_ArbitraryGraphLaplacian`](@ref) documentation for
more details on the `{-1,0,1}`-spectral properties of this Laplacian type.

See also the documentation for parent type [`_TypedLaplacian`](@ref) if seeking
justification for typing the `matrix` field as `Matrix{Int}` rather than
`AbstractMatrix{<:Integer}`.
"""
struct _ArbitraryGraphLaplacian <: _TypedLaplacian
    matrix::Matrix{Int}
end

"""
    _cast_to_typed_laplacian(L)

Classify the Laplacian matrix `L` and wrap it in a [`_TypedLaplacian`](@ref) object.

It is first verified that `L` is indeed a Laplacian matrix by
[`_assert_matrix_is_undirected_laplacian`](@ref), which throws a `DomainError` otherwise. It
is then classified based on any properties which may be exploited in computing data on its
`{-1,0,1}`-spectrum.

# Arguments
- `L::AbstractMatrix{<:Integer}`: the Laplacian matrix to classify.

# Returns
- `TL::_TypedLaplacian`: the Laplacian wrapped in a concrete
    [`_TypedLaplacian`](@ref) subtype associated with the category of the
    graph represented by `L`. In the case of complete graphs, `TL` also contains data on the
    (necessarily uniform) weight of each edge.

# Examples
Correctly recognizes the Laplacian matrix of the null graph:
```jldoctest; setup = :(using SDiagonalizability, Graphs)
julia> G = Graph(0)
{0, 0} undirected simple Int64 graph

julia> L = laplacian_matrix(G)
0×0 SparseArrays.SparseMatrixCSC{Int64, Int64} with 0 stored entries

julia> TL = SDiagonalizability._cast_to_typed_laplacian(L)
SDiagonalizability._NullGraphLaplacian(Matrix{Int64}(undef, 0, 0))

julia> isa(TL, SDiagonalizability._NullGraphLaplacian)
true
```

Correctly recognizes the Laplacian matrix of the empty graph on `3` nodes:
```jldoctest; setup = :(using SDiagonalizability, Graphs)
julia> G = Graph(3)
{3, 0} undirected simple Int64 graph

julia> L = laplacian_matrix(G)
3×3 SparseArrays.SparseMatrixCSC{Int64, Int64} with 0 stored entries:
 ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅

julia> TL = SDiagonalizability._cast_to_typed_laplacian(L)
SDiagonalizability._EmptyGraphLaplacian([0 0 0; 0 0 0; 0 0 0])

julia> isa(TL, SDiagonalizability._EmptyGraphLaplacian)
true
```

Correctly recognizes the Laplacian matrix of the complete graph on `4` nodes:
```jldoctest; setup = :(using SDiagonalizability, Graphs)
julia> G = complete_graph(4)
{4, 6} undirected simple Int64 graph

julia> L = laplacian_matrix(G)
4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 16 stored entries:
  3  -1  -1  -1
 -1   3  -1  -1
 -1  -1   3  -1
 -1  -1  -1   3

julia> TL = SDiagonalizability._cast_to_typed_laplacian(L)
SDiagonalizability._CompleteGraphLaplacian([3 -1 -1 -1; -1 3 -1 -1; -1 -1 3 -1; -1 -1 -1 3], -1)

julia> isa(TL, SDiagonalizability._CompleteGraphLaplacian)
true
```

Correctly recognizes the Laplacian matrix of this random generic graph:
```jldoctest; setup = :(using SDiagonalizability, Graphs)
julia> G = erdos_renyi(5, 8; seed=87)
{5, 8} undirected simple Int64 graph

julia> L = laplacian_matrix(G)
5×5 SparseArrays.SparseMatrixCSC{Int64, Int64} with 21 stored entries:
  3   ⋅  -1  -1  -1
  ⋅   2  -1  -1   ⋅
 -1  -1   4  -1  -1
 -1  -1  -1   4  -1
 -1   ⋅  -1  -1   3

julia> TL = SDiagonalizability._cast_to_typed_laplacian(L)
SDiagonalizability._ArbitraryGraphLaplacian([3 0 … -1 -1; 0 2 … -1 0; … ; -1 -1 … 4 -1; -1 0 … -1 3])

julia> isa(TL, SDiagonalizability._ArbitraryGraphLaplacian)
true
```
# Notes
Right now, the possible types of Laplacian matrices (or, to be more precise, the graphs they
represent) are:
- [`_NullGraphLaplacian`](@ref): the (unique) graph with no nodes.
- [`_EmptyGraphLaplacian`](@ref): any graph with no edges on `n ≥ 1` nodes.
- [`_CompleteGraphLaplacian`](@ref): any graph on `n ≥ 2` nodes where every pair of nodes is
    connected by an edge and all edges possess the same weight.
- [`_ArbitraryGraphLaplacian`](@ref): any graph on `n ≥ 3` nodes with no particular
    strcuture of relevance in the context of investigating `{-1,0,1}`-spectra.
"""
function _cast_to_typed_laplacian(L::AbstractMatrix{<:Integer})
    #= Verify that `L` is symmetric (thus representing an undirected graph) and has zero row
    sums (thus being a Laplacian matrix). =#
    _assert_matrix_is_undirected_laplacian(L)

    L_copy = Matrix{Int}(L) # Avoid shared mutability and cast to `Matrix{Int}`
    n = size(L_copy, 1)

    if n == 0 # The graph has no nodes
        TL = _NullGraphLaplacian(L_copy)
    elseif iszero(L_copy) # The graph has no edges
        TL = _EmptyGraphLaplacian(L_copy)
    #! format: off
    #= We have already verified symmetry, so equality of all strictly upper-triangular
    elements is a sufficient condition for `L` to represent a complete graph. =#
    #! format: on
    elseif allequal(L_copy[i, j] for i in 1:(n - 1) for j in (i + 1):n)
        weight = L_copy[1, 2]
        TL = _CompleteGraphLaplacian(L_copy, weight)
    else # The graph has no particular structure of interest
        TL = _ArbitraryGraphLaplacian(L_copy)
    end

    return TL
end
