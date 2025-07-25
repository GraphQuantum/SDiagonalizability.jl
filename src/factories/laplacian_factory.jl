# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    ClassifiedLaplacian

An abstract type representing a classified Laplacian matrix of an undirected graph.

# Interface
Concrete subtypes of `_TypedLaplacian` *must* implement the following fields:
- `matrix::Matrix{Int}`: the Laplacian matrix of the graph.
"""
abstract type ClassifiedLaplacian end

"""
    NullGraphLaplacian <: ClassifiedLaplacian

[TODO: Write here]
"""
struct NullGraphLaplacian <: ClassifiedLaplacian
    matrix::Matrix{Int}
end

"""
    EmptyGraphLaplacian <: ClassifiedLaplacian

[TODO: Write here]
"""
struct EmptyGraphLaplacian <: ClassifiedLaplacian
    matrix::Matrix{Int}
end

"""
    CompleteGraphLaplacian <: ClassifiedLaplacian

[TODO: Write here]
"""
struct CompleteGraphLaplacian <: ClassifiedLaplacian
    matrix::Matrix{Int}
    weight::Int
end

"""
    ArbitraryGraphLaplacian <: ClassifiedLaplacian

[TODO: Write here]
"""
struct ArbitraryGraphLaplacian <: ClassifiedLaplacian
    matrix::Matrix{Int}
end

"""
    classify_laplacian(L)

Classify the Laplacian matrix `L` and wrap it in a [`ClassifiedLaplacian`](@ref) object.

It is first verified that `L` is indeed a Laplacian matrix by
[`_assert_matrix_is_undirected_laplacian`](@ref), which throws a `DomainError` otherwise. It
is then classified based on any properties which may be exploited in computing data on its
`{-1,0,1}`-spectrum.

# Arguments
- `L::AbstractMatrix{<:Integer}`: the Laplacian matrix to classify.

# Returns
- `CL::ClassifiedLaplacian`: the Laplacian wrapped in a concrete
    [`ClassifiedLaplacian`](@ref) subtype associated with the category of the
    graph represented by `L`. In the case of complete graphs, `CL` also contains data on the
    (necessarily uniform) weight of each edge.

# Examples
Correctly recognizes the Laplacian matrix of the null graph:
```jldoctest
julia> using Graphs

julia> G = Graph(0)
{0, 0} undirected simple Int64 graph

julia> L = laplacian_matrix(G)
0×0 SparseArrays.SparseMatrixCSC{Int64, Int64} with 0 stored entries

julia> SDiagonalizability.classify_laplacian(L)
SDiagonalizability.NullGraphLaplacian(Matrix{Int64}(undef, 0, 0))
```

Correctly recognizes the Laplacian matrix of the empty graph on ``3`` nodes:
```jldoctest
julia> using Graphs

julia> G = Graph(3)
{3, 0} undirected simple Int64 graph

julia> L = laplacian_matrix(G)
3×3 SparseArrays.SparseMatrixCSC{Int64, Int64} with 0 stored entries:
 ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅

julia> SDiagonalizability.classify_laplacian(L)
SDiagonalizability.EmptyGraphLaplacian([0 0 0; 0 0 0; 0 0 0])
```

Correctly recognizes the Laplacian matrix of the complete graph on ``4`` nodes:
```jldoctest
julia> using Graphs

julia> G = complete_graph(4)
{4, 6} undirected simple Int64 graph

julia> L = laplacian_matrix(G)
4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 16 stored entries:
  3  -1  -1  -1
 -1   3  -1  -1
 -1  -1   3  -1
 -1  -1  -1   3

julia> SDiagonalizability.classify_laplacian(L)
SDiagonalizability.CompleteGraphLaplacian([3 -1 -1 -1; -1 3 -1 -1; -1 -1 3 -1; -1 -1 -1 3], 1)
```

Correctly recognizes the Laplacian matrix of this random generic graph:
```jldoctest
julia> using Graphs

julia> G = erdos_renyi(5, 8; seed=87)
{5, 8} undirected simple Int64 graph

julia> L = laplacian_matrix(G)
5×5 SparseArrays.SparseMatrixCSC{Int64, Int64} with 21 stored entries:
  3   ⋅  -1  -1  -1
  ⋅   2  -1  -1   ⋅
 -1  -1   4  -1  -1
 -1  -1  -1   4  -1
 -1   ⋅  -1  -1   3

julia> SDiagonalizability.classify_laplacian(L)
SDiagonalizability.ArbitraryGraphLaplacian([3 0 … -1 -1; 0 2 … -1 0; … ; -1 -1 … 4 -1; -1 0 … -1 3])
```
# Notes
Right now, the possible types of Laplacian matrices (or, to be more precise, the graphs they
represent) are:
- [`NullGraphLaplacian`](@ref): the (unique) graph with no nodes.
- [`EmptyGraphLaplacian`](@ref): any graph with no edges on ``n ≥ 1`` nodes.
- [`CompleteGraphLaplacian`](@ref): any graph on ``n ≥ 2`` nodes where every pair of nodes is
    connected by an edge and all edges possess the same weight.
- [`ArbitraryGraphLaplacian`](@ref): any graph on ``n ≥ 3`` nodes with no particular
    strcuture of relevance in the context of investigating ``{-1,0,1}``-spectra.
"""
function classify_laplacian(L::AbstractMatrix{<:Integer})
    #= Verify that `L` is symmetric (thus representing an undirected graph) and has zero row
    sums (thus being a Laplacian matrix). =#
    _assert_matrix_is_undirected_laplacian(L)

    L_copy = Matrix{Int}(L) # Avoid shared mutability and cast to `Matrix{Int}`
    n = size(L_copy, 1)

    if n == 0 # The graph has no nodes
        CL = NullGraphLaplacian(L_copy)
    elseif iszero(L_copy) # The graph has no edges
        CL = EmptyGraphLaplacian(L_copy)
    #! format: off
    #= We have already verified symmetry, so equality of all strictly upper-triangular
    elements is a sufficient condition for `L` to represent a complete graph. =#
    #! format: on
    elseif allequal(L_copy[i, j] for i in 1:(n - 1) for j in (i + 1):n)
        weight = -L_copy[1, 2]
        CL = CompleteGraphLaplacian(L_copy, weight)
    else # The graph has no particular structure of interest
        CL = ArbitraryGraphLaplacian(L_copy)
    end

    return CL
end
