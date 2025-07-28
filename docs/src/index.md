```@meta
CurrentModule = SDiagonalizability
```

```@raw html
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
```

## Overview

*SDiagonalizability.jl* implements the first non-naïve deterministic algorithm to minimize the ``S``-bandwidth of an undirected graph with integer edge weights.

Given an undirected, possibly weighted graph ``G`` and finite set of integers ``S \subset \mathbb{Z}``, ``G`` is said to be "``S``-diagonalizable" if there exists some diagonal matrix ``D`` and matrix ``P`` with all entries from ``S`` such that ``G``'s Laplacian matrix ``L(G) = PDP^{-1}``. If ``G`` is ``S``-diagonalizable, then its "``S``-bandwidth" is the minimum integer ``k \in \{1, 2, \ldots, |V(G)|\}`` such that there exists some diagonal matrix ``D`` and matrix ``P`` with all entries from ``S`` such that ``L(G) = PDP^{-1}`` and ``[P^\mathsf{T}P]_{i,j} = 0`` whenever ``|i - j| \geq k``; otherwise, its ``S``-bandwidth is simply ``\infty``.

For specific choices of ``S`` (namely ``\{-1, 1\}`` and ``\{-1, 0, 1\}``), the ``S``-bandwidth of a quantum network has been shown to be an indicator of high state transfer fidelity due to automorphic properties of the graph. As such, the nascent study of ``S``-diagonalizability and ``S``-bandwidth is of interest in the broader context of quantum information theory.

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

[TODO: Write here]

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

## Index

```@index
```
