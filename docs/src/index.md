```@meta
CurrentModule = SDiagonalizability
```

# SDiagonalizability.jl

*SDiagonalizability.jl* packages the first non-naïve deterministic algorithm to minimize the ``S``-bandwidth of a quantum network (or, to be more precise, its graph representation), written in Julia. Given some finite ``S \in \mathbb{Z}^n``, the ``S``-bandwidth of an (undirected) graph ``G`` with Laplacian ``L(G) \in \mathbb{R}^n`` is the minimum ``k \ge 1`` such that ``L(G) = PDP^{-1}`` for some diagonal ``D \in \mathbb{R}^{n \times n}`` and ``P \in S^{n \times n}`` with ``P^T P_{ij} = 0`` whenever ``|i - j| \ge k``. For specific choices of ``S`` (namely ``S = \{-1,1\}`` and ``S = \{-1,0,1\}``), a quantum network's ``S``-bandwidth has been shown to be an indicator of high state transfer fidelity due to automorphic properties of its graph representation.

!!! important "CURRENTLY UNDER DEVELOPMENT"
    We aim to release the initial stable version of *SDiagonalizability.jl* in mid-June 2025. The current version is a work-in-progress, with the main `SDiagonalizability` module not yet functional and comprehensive documentation/tests not yet available.

## Credits

Created by [Luis M. B. Varona](https://github.com/Luis-Varona), [Dr. Nathaniel Johnston](https://github.com/nathanieljohnston), and Dr. Sarah Plosker.

Insights from [Benjamin Talbot](https://github.com/Benjamin-Talbot), [Logan Pipes](https://github.com/logan-pipes), and [Dr. Liam Keliher](https://github.com/lkeliher) are gratefully acknowledged.

## Index

```@index
```
