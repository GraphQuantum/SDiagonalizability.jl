# SDiagonalizability.jl

<table>
  <tr>
    <td>
      <table>
        <tr>
          <td>Metadata</td>
          <td>
            <img src="https://img.shields.io/badge/version-v0.1.0--dev-slateblue" alt="Version">
            <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-darkorchid" alt="License: MIT"></a>
            <a href="https://github.com/JuliaDiff/BlueStyle"><img src="https://img.shields.io/badge/code%20style-blue-4495d1.svg" alt="Code Style: Blue"></a>
          </td>
        </tr>
        <tr>
          <td>Documentation</td>
          <td>
            <a href="https://graphquantum.github.io/SDiagonalizability.jl/stable"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Documentation of latest stable version"></a>
            <a href="https://graphquantum.github.io/SDiagonalizability.jl/dev"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Documentation of dev version"></a>
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
            <a href="https://codecov.io/gh/GraphQuantum/SDiagonalizability.jl"><img src="https://img.shields.io/codecov/c/gh/GraphQuantum/SDiagonalizability.jl?label=codecov" alt="Test coverage from codecov"></a>
          </td>
        </tr>
        <tr>
          <td>Static analysis with</td>
          <td>
            <a href="https://github.com/JuliaTesting/Aqua.jl"><img src="https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg" alt="Aqua QA"></a>
            <a href="https://github.com/aviatesk/JET.jl"><img src="https://img.shields.io/badge/%E2%9C%88%20tested%20with-JET.jl%EF%B8%8F-9cf" alt="JET static analysis"></a>
          </td>
        </tr>
      </table>
    </td>
    <td>
      <img src="https://github.com/GraphQuantum/SDiagonalizability.jl/raw/main/docs/src/assets/logo.jpg" alt="GraphQuantum logo by Madiha Waqar" width="170"><br>
      <sup><sub><em>GraphQuantum</em> logo by <a href="https://github.com/madihaahmed1">Madiha Waqar</a></sub></sup>
    </td>
  </tr>
</table>

*SDiagonalizability.jl* packages the first non-naïve deterministic algorithm to minimize the *S*-bandwidth of a quantum network (or, to be more precise, its graph representation), written in Julia. Given some finite *S* &isin; **Z**<sup>*n*</sup>, the *S*-bandwidth of an (undirected) graph *G* with Laplacian *L(G)* &isin; **R**<sup>*n*&times;*n*</sup> is the minimum *k* &ge; 1 such that *L(G)* = *PDP*<sup>-1</sup> for some diagonal *D* &isin; **R**<sup>*n*&times;*n*</sup> and *P* &isin; **S**<sup>*n*&times;*n*</sup> with [*P<sup>T</sup>P*]<sub>*ij*</sub> = 0 whenever |*i* - *j*| &ge; *k*. For specific choices of *S* (namely *S* = {-1,1} and *S* = {-1,0,1}), a quantum network's *S*-bandwidth has been shown to be an indicator of high state transfer fidelity due to automorphic properties of its graph representation.

**CURRENTLY UNDER DEVELOPMENT:** We aim to release the initial stable version of *SDiagonalizability.jl* in mid-June 2025. The current version is a work-in-progress, with the main `SDiagonalizability` module not yet functional and comprehensive documentation/tests not yet available.

## Credits

Created by [Luis M. B. Varona](https://github.com/Luis-Varona), [Dr. Nathaniel Johnston](https://github.com/nathanieljohnston), and Dr. Sarah Plosker.

Insights from [Benjamin Talbot](https://github.com/Benjamin-Talbot), [Logan Pipes](https://github.com/logan-pipes), and [Dr. Liam Keliher](https://github.com/lkeliher) are gratefully acknowledged.
