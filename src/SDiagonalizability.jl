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
using MolecularGraph: subgraph_monomorphisms

using PrecompileTools: @setup_workload, @compile_workload

include("utils.jl")
include("factories/laplacian_types.jl")
include("factories/orthogonality_properties.jl")

include("eigenvector_generation.jl")
include("laplacian_spectra.jl")
include("basis_search.jl")
include("s_bandwidth.jl")

export SpectrumIntegralResult, check_spectrum_integrality
# TODO: Add exports from `s_bandwidth.jl`

end
