# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SDiagonalizability

A dynamic algorithm to minimize or recognize the ``S``-bandwidth of an undirected graph.

Given an undirected, possibly weighted graph ``G`` and finite set of integers ``S`` ⊂ ``ℤ``,
``G`` is said to be "``S``-diagonalizable" if there exists some diagonal matrix ``D`` and
matrix ``P`` with all entries from ``S`` such that ``G``'s Laplacian matrix
``L(G) = PDP⁻¹``. If ``G`` is ``S``-diagonalizable, then its "``S``-bandwidth" is the
minimum integer ``k ∈ \\{1, 2, …, |V(G)|\\}`` such that there exists some diagonal matrix
``D`` and matrix ``P`` with all entries from ``S`` such that ``L(G) = PDP⁻¹`` and
``[PᵀP]ᵢⱼ = 0`` whenever ``|i - j| ≥ k``; otherwise, its ``S``-bandwidth is simply ``∞``.

For specific choices of ``S`` (namely ``\\{-1, 1\\}`` and ``\\{-1, 0, 1\\}``), the
``S``-bandwidth of a quantum network has been shown to be an indicator of high state
transfer fidelity due to automorphic properties of the graph. As such, the nascent study of
``S``-diagonalizability and ``S``-bandwidth is of interest in the broader context of quantum
information theory.

This package, therefore, implements the first algorithm beyond mere brute force to minimize
the ``S``-bandwidth of an undirected graph with integer edge weights. Capabilities also
exist for determining whether a graph has ``S``-bandwidth less than or equal to a fixed
integer ``k ≥ 1`` without necessarily caring about the true minimum value.

The full documentation is available at
[GitHub Pages](https://graphquantum.github.io/SDiagonalizability.jl/).
"""
module SDiagonalizability

using Graphs
using LinearAlgebra
using MatrixBandwidth

using Combinatorics: combinations, multiset_permutations
using DataStructures: OrderedDict, counter
using ElasticArrays: ElasticMatrix
using PrecompileTools: @setup_workload, @compile_workload

include("utils.jl")
include("types.jl")

include("factories/laplacian_factory.jl")
include("factories/orthogonality_factory.jl")

include("eigenvector_generation.jl")
include("laplacian_s_spectra.jl")
include("basis_search.jl")

include("core.jl")

include("startup.jl")

export s_bandwidth, has_s_bandwidth_at_most_k, is_s_diagonalizable

end
