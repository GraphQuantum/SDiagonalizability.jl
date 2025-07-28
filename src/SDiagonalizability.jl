# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SDiagonalizability

A dynamic algorithm to minimize the S-bandwidth of an undirected graph.

[TODO: Write here]
"""
module SDiagonalizability

using Combinatorics
using DataStructures
using ElasticArrays
using Graphs
using LinearAlgebra

include("utils.jl")
include("types.jl")

include("factories/laplacian_factory.jl")
# include("factories/orthogonality_factory.jl")

include("eigenvector_generation.jl")
include("laplacian_s_spectra.jl")
# include("basis_search.jl")

# include("core.jl")

# export minimize_s_bandwidths, has_s_bands_at_most_k, is_s_diagonalizable

end
