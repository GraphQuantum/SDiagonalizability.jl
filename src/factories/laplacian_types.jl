# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

# TODO: Add comments and docstrings to this file

abstract type _TypedLaplacian <: AbstractMatrix{Int} end

Base.size(L::_TypedLaplacian) = size(L.matrix)
Base.getindex(L::_TypedLaplacian, i::Int, j::Int) = L.matrix[i, j]
Base.eltype(L::_TypedLaplacian) = eltype(L.matrix)

struct _NullGraphLaplacian <: _TypedLaplacian
    matrix::Matrix{Int}
end

struct _EmptyGraphLaplacian <: _TypedLaplacian
    matrix::Matrix{Int}
end

struct _CompleteGraphLaplacian <: _TypedLaplacian
    matrix::Matrix{Int}
    weight::Int
end

struct _ArbitraryGraphLaplacian <: _TypedLaplacian
    matrix::Matrix{Int}
end

function _cast_to_typed_laplacian(L::AbstractMatrix{<:Integer})
    _assert_matrix_is_undirected_laplacian(L)
    L_copy = Matrix{Int}(L) # Avoid shared mutability and cast to `Matrix{Int}`
    n = size(L_copy, 1)

    if n == 0
        TL = _NullGraphLaplacian(L_copy)
    elseif iszero(L_copy)
        TL = _EmptyGraphLaplacian(L_copy)
    elseif allequal(L_copy[i, j] for i in 1:n, j in 1:n if i != j)
        weight = L_copy[1, 2]
        TL = _CompleteGraphLaplacian(L_copy, weight)
    else
        TL = _ArbitraryGraphLaplacian(L_copy)
    end

    return TL
end
