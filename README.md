<img src="https://github.com/GraphQuantum/SDiagonalizability.jl/raw/main/docs/src/assets/logo.jpg" alt="GraphQuantum logo by Madiha Waqar" height="175"/>
<sup><sub><em>GraphQuantum</em> logo by <a href="https://github.com/madihaahmed1">Madiha Waqar</a></sub></sup>

# SDiagonalizability

<table>
  <tr>
    <td>Metadata</td>
    <td>
      <img src="https://img.shields.io/badge/version-v0.1.0--dev-pink.svg" alt="Version">
      <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-A31F34.svg" alt="License: MIT"></a>
      <a href="https://github.com/JuliaDiff/BlueStyle"><img src="https://img.shields.io/badge/code%20style-blue-4495d1.svg" alt="Code Style: Blue"></a>
    </td>
  </tr>
  <tr>
    <td>Documentation</td>
    <td>
      <a href="https://graphquantum.github.io/SDiagonalizability.jl/stable/"><img src="https://img.shields.io/badge/docs-stable-darkgreen.svg" alt="Documentation of latest stable version"></a>
      <a href="https://graphquantum.github.io/SDiagonalizability.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-rebeccapurple.svg" alt="Documentation of dev version"></a>
    </td>
  </tr>
  <tr>
    <td>Continuous integration</td>
    <td>
      <a href="https://github.com/GraphQuantum/SDiagonalizability.jl/actions?query=workflow%3ACI+branch%3Amain"><img src="https://github.com/GraphQuantum/SDiagonalizability.jl/actions/workflows/CI.yml/badge.svg" alt="GitHub Workflow Status"></a>
    </td>
  </tr>
  <tr>
    <td>Code coverage</td>
    <td>
      <a href="https://codecov.io/gh/GraphQuantum/SDiagonalizability.jl"><img src="https://codecov.io/gh/GraphQuantum/SDiagonalizability.jl/branch/main/graph/badge.svg" alt="Test coverage from codecov"></a>
    </td>
    </tr>
    <tr>
      <td>Static analysis with</td>
      <td>
        <a href="https://github.com/JuliaTesting/Aqua.jl"><img src="https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg" alt="Aqua QA"></a>
        <a href="https://github.com/aviatesk/JET.jl"><img src="https://img.shields.io/badge/%E2%9C%88%20tested%20with-JET.jl%EF%B8%8F-9cf.svg" alt="JET static analysis"></a>
      </td>
    </tr>
</table>

## Overview

*SDiagonalizability.jl* implements the first non-naïve deterministic algorithm to minimize the *S*-bandwidth of an undirected graph with integer edge weights. Capabilities also exist for determining whether a graph has *S*-bandwidth less than or equal to a fixed integer *k* &ge; 1 without necessarily caring about the true minimum value.

Given an undirected, possibly weighted graph *G* and finite set of integers *S* &subset; **Z**, *G* is said to be "*S*-diagonalizable" if there exists some diagonal matrix *D* and matrix *P* with all entries from *S* such that *G*'s Laplacian matrix *L*(*G*) = *PDP*<sup>-1</sup>. If *G* is *S*-diagonalizable, then its "*S*-bandwidth" is the minimum integer *k* &isin; {1, 2, &hellip;, |*V*(*G*)|} such that there exists some diagonal matrix *D* and matrix *P* with all entries from *S* such that *L*(*G*) = *PDP*<sup>-1</sup> and [*P*<sup>T</sup>*P*]<sub>*i,j*</sub> = 0 whenever |*i* - *j*| &ge; *k*; otherwise, its *S*-bandwidth is simply &infin;.

For specific choices of *S* (namely {-1, 1} and {-1, 0, 1}), the *S*-bandwidth of a quantum network has been shown to be an indicator of high state transfer fidelity due to automorphic properties of the graph. As such, the nascent study of *S*-diagonalizability and *S*-bandwidth is of interest in the broader context of quantum information theory.

## Installation

The only prerequisite is a working Julia installation (v1.10 or later). First, enter Pkg mode by typing `]` in the Julia REPL, then run the following command:

```julia-repl
pkg> add https://github.com/GraphQuantum/SDiagonalizability.jl
```

Once *SDiagonalizability.jl* is added to Julia's General package registry, you will be able to install it more easily with:

```julia-repl
pkg> add SDiagonalizability
```

## Basic use

*SDiagonalizability.jl* offers capabilities for *S*-bandwidth minimization, *S*-bandwidth recognition, and plain old *S*-diagonalizability checking. To compute the *S*-bandwidth of a graph, you can use the `s_bandwidth` function:

```julia-repl
julia> using Graphs

julia> g = complete_graph(14)
{14, 91} undirected simple Int64 graph

julia> res_01neg = s_bandwidth(g, (-1, 0, 1))
Results of S-Bandwidth Minimization
 * S: (-1, 0, 1)
 * S-Bandwidth: 2
 * Graph Order: 14

julia> res_01neg.s_diagonalization
LinearAlgebra.Eigen{Int64, Int64, Matrix{Int64}, Vector{Int64}}
values:
14-element Vector{Int64}:
  0
 14
 14
 14
 14
 14
 14
 14
 14
 14
 14
 14
 14
 14
vectors:
14×14 Matrix{Int64}:
 1   0   0   0   0   1   1   1   1   1   1   0   0   0
 1   0   0   0   0  -1   0  -1   1   1   1   0   0   0
 1   0   0   0   0   0  -1   1   0   1   1   0   0   0
 1   0   0   0   0   0   0  -1   0   1  -1   0   0   0
 1   0   0   0   0   0   0   0  -1  -1   1   1   1   1
 1   0   0   0   0   0   0   0  -1  -1   1  -1  -1  -1
 1   1   1   0   0   0   0   0   0  -1   0   0  -1   1
 1  -1  -1   0   0   0   0   0   0  -1   0   0   0  -1
 1   0  -1   1   1   0   0   0   0   0  -1   0   0   0
 1   0   0  -1  -1   0   0   0   0   0  -1   0   0   0
 1   0   0   0  -1   0   0   0   0   0  -1   0   0   0
 1   0   1   0   1   0   0   0   0   0  -1   0   0   0
 1  -1  -1   0   0   0   0   0   0   0   0   0   0   1
 1   1   1   0   0   0   0   0   0   0   0   0   1  -1

julia> res_1neg = s_bandwidth(g, (-1, 1))
Results of S-Bandwidth Minimization
 * S: (-1, 1)
 * S-Bandwidth: 13
 * Graph Order: 14

julia> res_1neg.s_diagonalization
LinearAlgebra.Eigen{Int64, Int64, Matrix{Int64}, Vector{Int64}}
values:
14-element Vector{Int64}:
  0
 14
 14
 14
 14
 14
 14
 14
 14
 14
 14
 14
 14
 14
vectors:
14×14 Matrix{Int64}:
 1   1   1   1   1   1   1   1   1   1   1   1   1   1
 1  -1  -1  -1  -1  -1  -1   1   1   1   1   1   1   1
 1  -1  -1  -1   1   1   1  -1  -1  -1   1   1   1   1
 1  -1  -1  -1   1   1   1   1   1   1  -1  -1  -1   1
 1  -1  -1  -1   1   1   1   1   1   1   1   1   1  -1
 1  -1   1   1  -1  -1   1  -1  -1   1  -1  -1   1  -1
 1  -1   1   1  -1   1  -1  -1   1  -1  -1   1  -1  -1
 1  -1   1   1   1  -1  -1   1  -1  -1   1  -1  -1   1
 1   1  -1   1  -1  -1   1  -1   1  -1   1  -1  -1   1
 1   1  -1   1  -1   1  -1   1  -1  -1  -1  -1   1  -1
 1   1  -1   1   1  -1  -1  -1  -1   1  -1   1  -1  -1
 1   1   1  -1  -1  -1   1   1  -1  -1  -1   1  -1  -1
 1   1   1  -1  -1   1  -1  -1  -1   1   1  -1  -1   1
 1   1   1  -1   1  -1  -1  -1   1  -1  -1  -1   1  -1
```

Alternatively, to determine whether a graph has *S*-bandwidth less than or equal to some fixed integer *k* &ge; 1 without
necessarily caring about the true (minimum) value, you can take advantage of `has_s_bandwidth_at_most_k`:

```julia-repl
julia> using Graphs

julia> g = PetersenGraph()
{10, 15} undirected simple Int64 graph

julia> res_01neg = has_s_bandwidth_at_most_k(g, (-1, 0, 1), 4)
Results of S-Bandwidth Recognition
 * S: (-1, 0, 1)
 * S-Bandwidth Threshold k: 4
 * Has S-Bandwidth ≤ k: true
 * Graph Order: 10

julia> res_01neg.s_diagonalization
LinearAlgebra.Eigen{Int64, Int64, Matrix{Int64}, Vector{Int64}}
values:
10-element Vector{Int64}:
 0
 5
 5
 5
 5
 2
 2
 2
 2
 2
vectors:
10×10 Matrix{Int64}:
 1   1   1   0   0   1   1   1   1   1
 1  -1  -1   1   1   0   0   0   0   1
 1   0   1  -1  -1   0  -1   0  -1  -1
 1   0   0   0   1   0   0   0  -1  -1
 1   0  -1   0  -1   1   1   0   0   0
 1  -1   0  -1   0   0   0   1   1   0
 1   1   0  -1  -1  -1   0  -1   0   1
 1   1  -1   1   0   0  -1   0   0  -1
 1   0   0   1   0  -1   0   0   0   0
 1  -1   1   0   1   0   0  -1   0   0

julia> res_1neg = has_s_bandwidth_at_most_k(g, (-1, 1), 10)
Results of S-Bandwidth Recognition
 * S: (-1, 1)
 * S-Bandwidth Threshold k: 10
 * Has S-Bandwidth ≤ k: false
 * Graph Order: 10

julia> isnothing(res_1neg.s_diagonalization)
true
```

Lastly, the `is_s_diagonalizable` function can be used to simply determine whether a graph is *S*-diagonalizable:

```julia-repl
julia> using Graphs

julia> L = laplacian_matrix(complete_multipartite_graph([1, 1, 2, 2, 3]))
9×9 SparseArrays.SparseMatrixCSC{Int64, Int64} with 71 stored entries:
  8  -1  -1  -1  -1  -1  -1  -1  -1
 -1   8  -1  -1  -1  -1  -1  -1  -1
 -1  -1   7   ⋅  -1  -1  -1  -1  -1
 -1  -1   ⋅   7  -1  -1  -1  -1  -1
 -1  -1  -1  -1   7   ⋅  -1  -1  -1
 -1  -1  -1  -1   ⋅   7  -1  -1  -1
 -1  -1  -1  -1  -1  -1   6   ⋅   ⋅
 -1  -1  -1  -1  -1  -1   ⋅   6   ⋅
 -1  -1  -1  -1  -1  -1   ⋅   ⋅   6

julia> res_01neg = is_s_diagonalizable(L, (-1, 0, 1))
Results of S-Diagonalizability Check
 * S: (-1, 0, 1)
 * S-Diagonalizable: true
 * Graph Order: 9

julia> res_01neg.s_diagonalization
LinearAlgebra.Eigen{Int64, Int64, Matrix{Int64}, Vector{Int64}}
values:
9-element Vector{Int64}:
 0
 6
 6
 7
 7
 9
 9
 9
 9
vectors:
9×9 Matrix{Int64}:
 1   0   0   0   0   1   1   1   1
 1   0   0   0   0   0  -1  -1   1
 1   0   0   1   1  -1  -1   1  -1
 1   0   0  -1  -1  -1  -1   1  -1
 1   0   0  -1   1  -1   1  -1   0
 1   0   0   1  -1  -1   1  -1   0
 1   1   1   0   0   1   0   0   0
 1  -1   0   0   0   1   0   0   0
 1   0  -1   0   0   1   0   0   0

julia> res_1neg = is_s_diagonalizable(L, (-1, 1))
Results of S-Diagonalizability Check
 * S: (-1, 1)
 * S-Diagonalizable: false
 * Graph Order: 9

julia> isnothing(res_1neg.s_diagonalization)
true
```

(As demonstrated in this last example, you can supply the Laplacian matrix of a graph as an alternative to the `Graph` object itself from the *Graphs.jl* package.)

## Documentation

The full documentation is available at [GitHub Pages](https://graphquantum.github.io/SDiagonalizability.jl/). Documentation for methods and types is also available via the Julia REPL. [TODO: Write here]

## Citing

We encourage you to cite our work if you have used our algorithm in your research. Starring the *SDiagonalizability.jl* repository on GitHub is also appreciated.

The latest citation information may be found in the [CITATION.bib](https://raw.githubusercontent.com/GraphQuantum/SDiagonalizability.jl/main/CITATION.bib) file within the repository.

## Project status

We aim to release the first stable version of *SDiagonalizability.jl* sometime in early August 2025. The current codebase is a work-in-progress, with the main `SDiagonalizability` module not yet functional and comprehensive documentation/tests not yet available.

## Credits

Created by [Luis M. B. Varona](https://github.com/Luis-Varona), [Dr. Nathaniel Johnston](https://github.com/nathanieljohnston), and Dr. Sarah Plosker.

Insights from [Benjamin Talbot](https://github.com/Benjamin-Talbot), [Logan Pipes](https://github.com/logan-pipes), and [Dr. Liam Keliher](https://github.com/lkeliher) are gratefully acknowledged.

Additional thanks to [Madiha Waqar](https://github.com/madihaahmed1) for the *GraphQuantum* logo and [Luc Campbell](https://github.com/Luc-Campbell) for code review.
