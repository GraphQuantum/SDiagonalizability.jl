# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    minimize_s_bandwidth(g::AbstractGraph, S)
    minimize_s_bandwidth(L::AbstractMatrix{<:Integer}, S)

[TODO: Write here]
"""
function minimize_s_bandwidth(g::AbstractGraph, S::Tuple{Vararg{Int}})
    _assert_graph_has_defined_s_bandwidth(g)
    return minimize_s_bandwidth(laplacian_matrix(g), S)
end

function minimize_s_bandwidth(L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Int}})
    _assert_matrix_is_undirected_laplacian(L)
    # TODO: Write here
    return nothing
end

"""
    has_s_bandwidth_at_most_k(g::AbstractGraph, S, k)
    has_s_bandwidth_at_most_k(L::AbstractMatrix{<:Integer}, S, k)

[TODO: Write here]
"""
function has_s_bandwidth_at_most_k(g::AbstractGraph, S::Tuple{Vararg{Int}}, k::Int)
    _assert_graph_has_defined_s_bandwidth(g)
    return has_s_bandwidth_at_most_k(laplacian_matrix(g), S, k)
end

function has_s_bandwidth_at_most_k(
    L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Int}}, k::Int
)
    _assert_matrix_is_undirected_laplacian(L)
    # TODO: Write here
    return nothing
end

"""
    is_s_diagonalizable(g::AbstractGraph, S)
    is_s_diagonalizable(L::AbstractMatrix{<:Integer}, S)

[TODO: Write here]
"""
function is_s_diagonalizable(g::AbstractGraph, S::Tuple{Vararg{Int}})
    _assert_graph_has_defined_s_bandwidth(g)
    return is_s_diagonalizable(laplacian_matrix(g), S)
end

function is_s_diagonalizable(L::AbstractMatrix{<:Integer}, S::Tuple{Vararg{Int}})
    _assert_matrix_is_undirected_laplacian(L)
    return has_s_bandwidth_at_most_k(L, S, size(L, 1))
end
