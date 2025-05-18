# SDiagonalizability.jl

![Version](https://img.shields.io/badge/version-v0.1.0--dev-slateblue)
![License: MIT](https://img.shields.io/badge/License-MIT-darkorchid)
[![Build Status](https://github.com/GraphQuantum/SDiagonalizability.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GraphQuantum/SDiagonalizability.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)

*SDiagonalizability.jl* packages the first non-naïve deterministic algorithm to minimize the *S*-bandwidth of a quantum network (or, to be more precise, its graph representation), written in Julia. Given some finite *S* &isin; **Z**<sup>*n*</sup>, the *S*-bandwidth of an (undirected) graph *G* with Laplacian *L(G)* &isin; **R**<sup>*n*&times;*n*</sup> is the minimum *k* &ge; 1 such that *L(G)* = *PDP*<sup>-1</sup> for some diagonal *D* &isin; **R**<sup>*n*&times;*n*</sup> and *P* &isin; **S**<sup>*n*&times;*n*</sup> with [*P<sup>T</sup>P*]<sub>*ij*</sub> = 0 whenever |*i* - *j*| &ge; *k*. For specific choices of *S* (namely *S* = {-1,1} and *S* = {-1,0,1}), a quantum network's *S*-bandwidth has been shown to be an indicator of high state transfer fidelity due to automorphic properties of its graph representation.

**(CURRENTLY UNDER DEVELOPMENT)**

## Credits

Created by [Luis M. B. Varona](https://github.com/Luis-Varona), [Dr. Nathaniel Johnston](https://github.com/nathanieljohnston), and Dr. Sarah Plosker.

Insights from [Benjamin Talbot](https://github.com/Benjamin-Talbot), [Logan Pipes](https://github.com/logan-pipes), and [Dr. Liam Keliher](https://github.com/lkeliher) are gratefully acknowledged.
