# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

module SDiagonalizability

using Graphs
using LinearAlgebra

using Combinatorics: combinations, multiset_permutations
using DataStructures: OrderedDict, counter
using ElasticArrays: ElasticMatrix
using MolecularGraph: nodesubgraph_is_isomorphic

using PrecompileTools: @setup_workload, @compile_workload

include("utils.jl")
include("laplacian_types.jl")
include("eigenvector_generation.jl")
include("laplacian_spectra.jl")
include("basis_search.jl")
include("s_bandwidth.jl")

# TODO: Add exports

end
