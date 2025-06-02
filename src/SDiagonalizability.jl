# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SDiagonalizability

A dynamic algorithm to test a graph for *S*-diagonalizability and minimize its *S*-bandwidth.

*SDiagonalizability.jl* packages the first non-naïve deterministic algorithm to minimize the
*S*-bandwidth of a quantum network (or, to be more precise, its graph representation),
written in Julia. Given some finite *S* ∈ **Z**ⁿ, the *S*-bandwidth of an (undirected) graph
*G* with Laplacian *L(G)* ∈ **R**ⁿˣⁿ is the minimum integer *k* ≥ 1 such that *L(G)* =
*PDP*⁻¹ for some diagonal *D* ∈ **R**ⁿˣⁿ and *P* ∈ *S*ⁿˣⁿ with [*P*ᵀ*P*]ᵢⱼ = 0 whenever
|*i* - *j*| ≥ *k*. For specific choices of *S* (namely *S* = {-1,1} and *S* = {-1,0,1}), a
quantum network's *S*-bandwidth has been shown to be an indicator of high state transfer
fidelity due to automorphic properties of its graph representation.

[Full documentation](https://graphquantum.github.io/SDiagonalizability.jl/dev/) is available
for the latest development version of this package.
"""
module SDiagonalizability

using Graphs
using LinearAlgebra

using Combinatorics: combinations, multiset_permutations
using DataStructures: OrderedDict, counter
using ElasticArrays: ElasticMatrix
using MolecularGraph: subgraph_monomorphisms

# TODO: Create and include `precompile_workload.jl`
using PrecompileTools: @setup_workload, @compile_workload

# TODO: Check if we really need the module preface for cross-referencing private functions

include("utils.jl")
include("factories/laplacian_types.jl")
include("factories/orthogonality_properties.jl")

include("eigenvector_generation.jl")
include("laplacian_spectra.jl")
include("basis_search.jl")
include("s_bandwidth.jl")

# Exports from `laplacian_spectra.jl`
export LaplacianSpectrum01Neg, laplacian_spectra_01neg
export SpectrumIntegralResult, check_spectrum_integrality

# Exports from `s_bandwidth.jl`
export s_bandwidth

end
