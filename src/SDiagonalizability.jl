# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

module SDiagonalizability

using Combinatorics: combinations, multiset_permutations
using DataStructures: OrderedDict, counter
using ElasticArrays: ElasticMatrix
using Graphs: SimpleGraph, is_bipartite, is_connected
using LinearAlgebra
using MolecularGraph: nodesubgraph_is_isomorphic
using PrecompileTools: @setup_workload, @compile_workload

include("error_messages.jl")
include("eigenvector_generators.jl")
include("laplacian_spectra.jl")
include("k_orthogonalizability.jl")
include("s_bandwidth.jl")

end
